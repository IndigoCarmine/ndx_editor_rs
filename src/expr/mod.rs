//! The selection expression language.
//!
//! ```text
//! or        = and , { "|" , and } ;
//! and       = unary , { ( "&" | "\" ) , unary } ;
//! unary     = "!" , unary
//!           | "bonded" , [ integer , "of" ] , unary
//!           | "within" , number , "of" , unary
//!           | primary ;
//! primary   = "(" , or , ")" | predicate | group_ref ;
//! predicate = ( "name" | "resname" | "resid" | "resnr" | "element" | "chain" | "type" ) , arg ;
//! group_ref = integer | ( quoted | bare ) , [ "#" , integer ] ;
//! ```
//!
//! `0 & !1` means what it means in `gmx make_ndx`; everything else is a superset.
//!
//! The structure-dependent parts (`name`, `resid`, `bonded`, `within`, ...) parse today but fail
//! at evaluation with a specific "needs a structure file" error, because the `.gro` / `.top`
//! readers are not written yet. Fixing the grammar now means adding them later changes no syntax.

pub mod eval;
pub mod lex;
pub mod parse;

use std::ops::Range;

use crate::glob::Pattern;
use crate::model::{GroupRef, IndexFile};

pub use crate::glob::Pattern as GlobPattern;
pub use eval::EvalCtx;
pub use parse::parse_expr;

pub type Span = Range<usize>;

#[derive(Clone, Debug, PartialEq)]
pub enum Expr {
    Ref(GroupRef),
    /// Needs a structure (or, for `type`, a topology).
    Pred { pred: Pred, span: Span },
    /// Atoms bonded to the operand, `depth` hops out. Needs a bond graph.
    Bonded {
        of: Box<Expr>,
        depth: u8,
        span: Span,
    },
    /// Atoms within `radius` nm of the operand. Needs coordinates.
    Within {
        radius: f32,
        of: Box<Expr>,
        span: Span,
    },
    Not(Box<Expr>),
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
    /// `A \ B`, and what `A & !B` is rewritten into.
    Diff(Box<Expr>, Box<Expr>),
}

/// A per-atom predicate. Every variant needs data an `.ndx` file does not carry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Pred {
    Name(Vec<Pattern>),
    ResName(Vec<Pattern>),
    ResId(Vec<(i32, i32)>),
    Element(Vec<Pattern>),
    Chain(Vec<char>),
    /// Force-field atom type. Needs a topology.
    Type(Vec<Pattern>),
    /// `[ moleculetype ]` name. Needs a topology — and it is the closest thing GROMACS has to a
    /// chain, since a `.gro` carries no chain column.
    Molecule(Vec<Pattern>),
    /// Literal atom numbers: `atomid 116`, `atomid 1-10,15`. The one predicate that needs nothing
    /// loaded — it is make_ndx's `a 1-10`, usable inside an expression.
    AtomId(crate::atomset::AtomSet),
}

impl Pred {
    pub fn keyword(&self) -> &'static str {
        match self {
            Pred::Name(_) => "name",
            Pred::ResName(_) => "resname",
            Pred::ResId(_) => "resid",
            Pred::Element(_) => "element",
            Pred::Chain(_) => "chain",
            Pred::Type(_) => "type",
            Pred::Molecule(_) => "molecule",
            Pred::AtomId(_) => "atomid",
        }
    }

    fn arg_str(&self) -> String {
        match self {
            Pred::Name(p)
            | Pred::ResName(p)
            | Pred::Element(p)
            | Pred::Type(p)
            | Pred::Molecule(p) => p.iter().map(Pattern::as_str).collect::<Vec<_>>().join("_"),
            Pred::ResId(rs) => rs
                .iter()
                .map(|(lo, hi)| {
                    if lo == hi {
                        lo.to_string()
                    } else {
                        format!("{lo}-{hi}")
                    }
                })
                .collect::<Vec<_>>()
                .join(","),
            Pred::Chain(cs) => cs.iter().collect(),
            Pred::AtomId(set) => crate::atomset::format_ranges(set),
        }
    }
}

/// Binding strength, for parenthesizing auto-generated names.
fn prec(e: &Expr) -> u8 {
    match e {
        Expr::Or(..) => 1,
        Expr::And(..) | Expr::Diff(..) => 2,
        Expr::Not(_) | Expr::Bonded { .. } | Expr::Within { .. } => 3,
        Expr::Ref(_) | Expr::Pred { .. } => 4,
    }
}

/// The make_ndx-style name for a result group: `0 & !1` -> `System_&_!Protein`.
pub fn auto_name(ndx: &IndexFile, e: &Expr) -> String {
    name_at(ndx, e, 0)
}

fn name_at(ndx: &IndexFile, e: &Expr, parent_prec: u8) -> String {
    let s = match e {
        // A glob is named after the pattern, not after whichever groups it happened to match.
        Expr::Ref(r) => match ndx.resolve_all(r) {
            Ok(ids) if ids.len() == 1 => sanitize(&ndx.groups[ids[0]].name),
            _ => sanitize(&r.to_string()),
        },
        Expr::Pred { pred, .. } => format!("{}_{}", pred.keyword(), sanitize(&pred.arg_str())),
        Expr::Bonded { of, depth, .. } => {
            let inner = name_at(ndx, of, 0);
            if *depth == 1 {
                format!("bonded({inner})")
            } else {
                format!("bonded{depth}({inner})")
            }
        }
        Expr::Within { radius, of, .. } => {
            let inner = name_at(ndx, of, 0);
            format!("within{radius}({inner})")
        }
        Expr::Not(x) => format!("!{}", name_at(ndx, x, prec(e))),
        Expr::And(a, b) => format!(
            "{}_&_{}",
            name_at(ndx, a, prec(e)),
            name_at(ndx, b, prec(e))
        ),
        Expr::Or(a, b) => format!(
            "{}_|_{}",
            name_at(ndx, a, prec(e)),
            name_at(ndx, b, prec(e))
        ),
        // Render a rewritten difference the way make_ndx would have shown it.
        Expr::Diff(a, b) => format!(
            "{}_&_!{}",
            name_at(ndx, a, prec(e)),
            name_at(ndx, b, 3)
        ),
    };
    if prec(e) < parent_prec {
        format!("({s})")
    } else {
        s
    }
}

fn sanitize(s: &str) -> String {
    s.chars()
        .map(|c| if c.is_whitespace() { '_' } else { c })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Group;

    fn ndx() -> IndexFile {
        IndexFile {
            groups: vec![
                Group::new("System", (1..=10).collect()),
                Group::new("Protein", vec![1, 2, 3]),
                Group::new("SOL", vec![4, 5]),
                Group::new("Water and ions", vec![4, 5, 6]),
            ],
        }
    }

    fn name_of(src: &str) -> String {
        let f = ndx();
        auto_name(&f, &parse_expr(src).unwrap())
    }

    #[test]
    fn make_ndx_style_names() {
        assert_eq!(name_of("0 & !1"), "System_&_!Protein");
        assert_eq!(name_of("1 | 2"), "Protein_|_SOL");
        assert_eq!(name_of("1 \\ 2"), "Protein_&_!SOL");
    }

    #[test]
    fn parenthesizes_when_precedence_demands() {
        assert_eq!(name_of("(0 | 1) & 2"), "(System_|_Protein)_&_SOL");
        assert_eq!(name_of("0 | 1 & 2"), "System_|_Protein_&_SOL");
    }

    #[test]
    fn whitespace_in_names_becomes_underscores() {
        assert_eq!(name_of("\"Water and ions\""), "Water_and_ions");
    }

    #[test]
    fn structure_expressions_get_names_too() {
        assert_eq!(name_of("element H & bonded 1"), "element_H_&_bonded(Protein)");
        assert_eq!(name_of("within 0.5 of 1"), "within0.5(Protein)");
    }

    #[test]
    fn glob_matching() {
        assert!(Pattern::new("CA").matches("CA"));
        assert!(!Pattern::new("CA").matches("CB"));
        assert!(Pattern::new("H*").matches("HW1"));
        assert!(Pattern::new("H*").matches("H"));
        assert!(!Pattern::new("H*").matches("OW"));
        assert!(Pattern::new("?W").matches("OW"));
        assert!(Pattern::new("*").matches("anything"));
        assert!(Pattern::new("C*A").matches("CBBA"));
        assert!(!Pattern::new("C*A").matches("CBBB"));
    }
}
