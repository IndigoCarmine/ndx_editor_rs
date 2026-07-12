use std::cell::{OnceCell, RefCell};
use std::collections::HashMap;
use std::rc::Rc;

use crate::atomset::AtomSet;
use crate::error::{NdxError, Result};
use crate::expr::{Expr, Pattern, Pred};
use crate::model::{AtomId, IndexFile};
use crate::system::SystemCtx;
use crate::universe::{self, Universe, UniverseSpec};

pub struct EvalCtx<'a> {
    pub ndx: &'a IndexFile,
    spec: &'a UniverseSpec,
    system: &'a SystemCtx,
    universe: OnceCell<Universe>,
    /// A group referenced five times is sorted once.
    cache: RefCell<HashMap<usize, Rc<AtomSet>>>,
    warnings: RefCell<Vec<String>>,
}

impl<'a> EvalCtx<'a> {
    pub fn new(ndx: &'a IndexFile, spec: &'a UniverseSpec, system: &'a SystemCtx) -> Self {
        EvalCtx {
            ndx,
            spec,
            system,
            universe: OnceCell::new(),
            cache: RefCell::new(HashMap::new()),
            warnings: RefCell::new(Vec::new()),
        }
    }

    /// Rewrite, then evaluate. Callers should use this rather than [`EvalCtx::eval`] directly.
    pub fn run(&self, e: &Expr) -> Result<Rc<AtomSet>> {
        let e = rewrite(e.clone(), self.rewrite_is_sound());
        self.eval(&e)
    }

    /// Diagnostics accumulated during evaluation, for the caller to print on stderr.
    pub fn take_warnings(&self) -> Vec<String> {
        std::mem::take(&mut self.warnings.borrow_mut())
    }

    /// `A & !B == A \ B` only holds when the universe contains every atom in the file. That is
    /// true when we derive it (`Auto`) or when `--natoms` was validated against the file, but a
    /// user-chosen `--universe GROUP` may deliberately be smaller.
    fn rewrite_is_sound(&self) -> bool {
        !matches!(self.spec, UniverseSpec::Group(_))
    }

    fn universe(&self) -> Result<&Universe> {
        // Resolved on first use, so an expression without a bare `!` never touches it — and
        // therefore never emits the "guessed universe" warning.
        if let Some(u) = self.universe.get() {
            return Ok(u);
        }
        let u = universe::resolve(self.ndx, self.spec, self.system)?;
        if let Some(w) = &u.warning {
            self.warnings.borrow_mut().push(w.clone());
        }
        Ok(self.universe.get_or_init(|| u))
    }

    fn group(&self, id: usize) -> Result<Rc<AtomSet>> {
        if let Some(s) = self.cache.borrow().get(&id) {
            return Ok(Rc::clone(s));
        }
        let set = Rc::new(self.ndx.get(id)?.to_set());
        self.cache.borrow_mut().insert(id, Rc::clone(&set));
        Ok(set)
    }

    pub fn eval(&self, e: &Expr) -> Result<Rc<AtomSet>> {
        Ok(match e {
            // A glob stands for the union of every group it matches, so `Fiber* | Alkyl*` reads
            // exactly the way it looks. An id or an exact name resolves to one group, and the
            // union of one set is itself.
            Expr::Ref(r) => {
                let ids = self.ndx.resolve_all(r)?;
                let mut it = ids.into_iter();
                let first = self.group(it.next().expect("resolve_all is never empty"))?;
                it.try_fold(first, |acc, id| {
                    let next = self.group(id)?;
                    Ok::<_, NdxError>(Rc::new(acc.union(&next)))
                })?
            }

            Expr::Not(x) => {
                let inner = self.eval(x)?;
                Rc::new(inner.complement(&self.universe()?.set))
            }

            Expr::And(a, b) => {
                let (a, b) = (self.eval(a)?, self.eval(b)?);
                Rc::new(a.intersection(&b))
            }

            Expr::Or(a, b) => {
                let (a, b) = (self.eval(a)?, self.eval(b)?);
                Rc::new(a.union(&b))
            }

            Expr::Diff(a, b) => {
                let (a, b) = (self.eval(a)?, self.eval(b)?);
                Rc::new(a.difference(&b))
            }

            // Everything below needs data an .ndx file does not carry. When the system was not
            // loaded, say exactly which flag would have supplied it.
            Expr::Pred { pred, span } => Rc::new(self.eval_pred(pred, span)?),

            Expr::Bonded { of, depth, span } => {
                let seed = self.eval(of)?;
                let graph = self.bonds().ok_or_else(|| missing(
                    "bonded",
                    "bond information",
                    "pass -p topol.top for real bonds, or -s conf.gro to estimate them \
                     from interatomic distances",
                    span,
                ))?;
                Rc::new(graph.expand(&seed, *depth))
            }

            Expr::Within { radius, of, span } => {
                let seed = self.eval(of)?;
                let s = self.structure().ok_or_else(|| {
                    missing("within", "atom coordinates", "pass -s conf.gro", span)
                })?;
                Rc::new(crate::spatial::within(s, &seed, *radius, self.system.pbc))
            }
        })
    }

    fn structure(&self) -> Option<&crate::structure::Structure> {
        self.system.structure.as_ref()
    }

    fn topology(&self) -> Option<&crate::topology::Topology> {
        self.system.topology.as_ref()
    }

    fn bonds(&self) -> Option<&crate::bonds::BondGraph> {
        self.system.bonds.as_ref()
    }

    /// Sweep a per-atom predicate over the whole system.
    fn eval_pred(&self, pred: &Pred, span: &crate::expr::Span) -> Result<AtomSet> {
        let want_structure = || {
            self.structure().ok_or_else(|| {
                missing(
                    pred.keyword(),
                    "a structure file",
                    "pass -s conf.gro",
                    span,
                )
            })
        };
        let want_topology = || {
            self.topology().ok_or_else(|| {
                missing(pred.keyword(), "a topology", "pass -p topol.top", span)
            })
        };

        // The one predicate that needs nothing loaded: it *is* the atom set.
        if let Pred::AtomId(set) = pred {
            return Ok(set.clone());
        }

        let hits: Vec<AtomId> = match pred {
            Pred::AtomId(_) => unreachable!("handled above"),

            Pred::Name(pats) => {
                let s = want_structure()?;
                sweep(s.natoms(), |a| {
                    s.name(a).is_some_and(|n| any_match(pats, n))
                })
            }

            Pred::ResName(pats) => {
                let s = want_structure()?;
                sweep(s.natoms(), |a| {
                    s.resname(a).is_some_and(|n| any_match(pats, n))
                })
            }

            Pred::Element(pats) => {
                let s = want_structure()?;
                // Element symbols are stored uppercase, so `element h` works too.
                let upper: Vec<Pattern> =
                    pats.iter().map(|p| Pattern::new(p.as_str().to_ascii_uppercase())).collect();
                sweep(s.natoms(), |a| {
                    s.element(a).is_some_and(|e| any_match(&upper, e))
                })
            }

            Pred::ResId(ranges) => {
                let s = want_structure()?;
                sweep(s.natoms(), |a| {
                    s.resid(a)
                        .is_some_and(|r| ranges.iter().any(|(lo, hi)| r >= *lo && r <= *hi))
                })
            }

            Pred::Type(pats) => {
                let t = want_topology()?;
                sweep(t.natoms, |a| {
                    t.atomtype(a).is_some_and(|ty| any_match(pats, ty))
                })
            }

            Pred::Molecule(pats) => {
                let t = want_topology()?;
                sweep(t.natoms, |a| {
                    t.molname(a).is_some_and(|m| any_match(pats, m))
                })
            }

            // A .gro has no chain column, and a .top has no chains at all — GROMACS splits systems
            // into molecules, not chains. Point at the thing that actually exists.
            Pred::Chain(_) => {
                return Err(NdxError::Other(
                    "`chain` is not available: a .gro file has no chain column\n\
                     hint: GROMACS splits a system into molecules, so `molecule <name>` (with \
                     -p topol.top) is usually what you want"
                        .into(),
                ));
            }
        };

        Ok(AtomSet::from_sorted_unique(hits))
    }
}

/// Every 1-based atom id satisfying `f`, in order — so the result is sorted by construction.
fn sweep(natoms: u32, mut f: impl FnMut(AtomId) -> bool) -> Vec<AtomId> {
    (1..=natoms).filter(|a| f(*a)).collect()
}

fn any_match(pats: &[Pattern], s: &str) -> bool {
    pats.iter().any(|p| p.matches(s))
}

fn missing(
    feature: &str,
    needs: &'static str,
    hint: &'static str,
    span: &crate::expr::Span,
) -> NdxError {
    NdxError::NeedsSystem {
        feature: feature.to_string(),
        needs,
        hint,
        span: span.clone(),
    }
}

/// Turn `A & !B` (and `!B & A`) into `A \ B`.
///
/// Set-identical whenever the universe contains A, and worth doing for two reasons: the
/// complement is never materialized, and the most common `!` expression by far stops triggering
/// the "I had to guess the universe" warning.
pub fn rewrite(e: Expr, enabled: bool) -> Expr {
    if !enabled {
        return e;
    }
    match e {
        Expr::And(a, b) => {
            let (a, b) = (rewrite(*a, true), rewrite(*b, true));
            match (a, b) {
                // Only when the *other* side is not itself a complement: `!A & !B` has no
                // universe-free form.
                (x, Expr::Not(y)) if !matches!(x, Expr::Not(_)) => {
                    Expr::Diff(Box::new(x), y)
                }
                (Expr::Not(y), x) if !matches!(x, Expr::Not(_)) => {
                    Expr::Diff(Box::new(x), y)
                }
                (a, b) => Expr::And(Box::new(a), Box::new(b)),
            }
        }
        Expr::Or(a, b) => Expr::Or(Box::new(rewrite(*a, true)), Box::new(rewrite(*b, true))),
        Expr::Diff(a, b) => Expr::Diff(Box::new(rewrite(*a, true)), Box::new(rewrite(*b, true))),
        Expr::Not(x) => Expr::Not(Box::new(rewrite(*x, true))),
        Expr::Bonded { of, depth, span } => Expr::Bonded {
            of: Box::new(rewrite(*of, true)),
            depth,
            span,
        },
        Expr::Within { radius, of, span } => Expr::Within {
            radius,
            of: Box::new(rewrite(*of, true)),
            span,
        },
        leaf @ (Expr::Ref(_) | Expr::Pred { .. }) => leaf,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expr::parse_expr;
    use crate::model::Group;

    fn ndx() -> IndexFile {
        IndexFile {
            groups: vec![
                Group::new("System", (1..=10).collect()),
                Group::new("Protein", vec![1, 2, 3]),
                Group::new("SOL", vec![3, 4, 5]),
            ],
        }
    }

    fn eval(src: &str) -> Result<Vec<u32>> {
        eval_with(src, UniverseSpec::Auto).map(|(v, _)| v)
    }

    fn eval_with(src: &str, spec: UniverseSpec) -> Result<(Vec<u32>, Vec<String>)> {
        let f = ndx();
        let system = SystemCtx::default();
        let cx = EvalCtx::new(&f, &spec, &system);
        let set = cx.run(&parse_expr(src).unwrap())?;
        Ok((set.as_slice().to_vec(), cx.take_warnings()))
    }

    #[test]
    fn intersection() {
        assert_eq!(eval("1 & 2").unwrap(), [3]);
    }

    #[test]
    fn union() {
        assert_eq!(eval("1 | 2").unwrap(), [1, 2, 3, 4, 5]);
    }

    #[test]
    fn difference_via_and_not() {
        assert_eq!(eval("1 & !2").unwrap(), [1, 2]);
    }

    #[test]
    fn difference_via_backslash() {
        assert_eq!(eval("1 \\ 2").unwrap(), [1, 2]);
    }

    #[test]
    fn complement_uses_the_system_group() {
        assert_eq!(eval("!1").unwrap(), [4, 5, 6, 7, 8, 9, 10]);
    }

    #[test]
    fn by_name_and_by_id_agree() {
        assert_eq!(eval("Protein & !SOL").unwrap(), eval("1 & !2").unwrap());
    }

    #[test]
    fn parens() {
        assert_eq!(eval("(1 | 2) & 0").unwrap(), [1, 2, 3, 4, 5]);
    }

    #[test]
    fn empty_result_is_not_an_error_here() {
        assert!(eval("1 & !1").unwrap().is_empty());
    }

    /// The point of the rewrite: `A & !B` must not drag the universe in, so a file without a
    /// `System` group produces no warning for it — while a bare `!B` still does.
    #[test]
    fn and_not_emits_no_universe_warning() {
        let f = IndexFile {
            groups: vec![
                Group::new("Protein", vec![1, 2, 3]),
                Group::new("SOL", vec![3, 4, 5]),
            ],
        };
        let spec = UniverseSpec::Auto;
        let system = SystemCtx::default();

        let cx = EvalCtx::new(&f, &spec, &system);
        cx.run(&parse_expr("0 & !1").unwrap()).unwrap();
        assert!(cx.take_warnings().is_empty(), "A & !B should not need a universe");

        let cx = EvalCtx::new(&f, &spec, &system);
        cx.run(&parse_expr("!1").unwrap()).unwrap();
        assert_eq!(cx.take_warnings().len(), 1, "a bare ! must warn");
    }

    /// With a deliberately undersized `--universe`, `A & !B` is *not* `A \ B`, so the rewrite
    /// must stay off and the real complement must be taken.
    #[test]
    fn rewrite_is_disabled_for_an_explicit_universe() {
        use crate::model::GroupRef;
        let spec = UniverseSpec::Group(GroupRef::name("SOL"));
        // Protein = {1,2,3}; universe = SOL = {3,4,5}; !Protein within it = {4,5}.
        // So Protein & !Protein-complement... concretely: 1 & !1 under universe SOL.
        let (atoms, _) = eval_with("1 & !1", spec).unwrap();
        assert!(atoms.is_empty());

        // 0 & !1 : System={1..10} ∩ (SOL \ Protein) = {4,5}, NOT System \ Protein = {4..10}.
        let spec = UniverseSpec::Group(GroupRef::name("SOL"));
        let (atoms, _) = eval_with("0 & !1", spec).unwrap();
        assert_eq!(atoms, [4, 5]);
    }

    #[test]
    fn natoms_widens_the_universe() {
        let (atoms, warnings) = eval_with("!0", UniverseSpec::Natoms(12)).unwrap();
        assert_eq!(atoms, [11, 12]);
        assert!(warnings.is_empty());
    }

    /// A glob stands for the union of the groups it matches.
    #[test]
    fn a_glob_unions_its_matches() {
        let f = IndexFile {
            groups: vec![
                Group::new("Fiber1", vec![1, 2]),
                Group::new("Fiber2", vec![3, 4]),
                Group::new("SOL", vec![5]),
            ],
        };
        let spec = UniverseSpec::Auto;
        let system = SystemCtx::default();
        let cx = EvalCtx::new(&f, &spec, &system);

        let set = cx.run(&parse_expr("Fiber*").unwrap()).unwrap();
        assert_eq!(set.as_slice(), [1, 2, 3, 4]);

        // And composes with the operators.
        let set = cx.run(&parse_expr("Fiber* | SOL").unwrap()).unwrap();
        assert_eq!(set.as_slice(), [1, 2, 3, 4, 5]);

        let set = cx.run(&parse_expr("Fiber* & !Fiber2").unwrap()).unwrap();
        assert_eq!(set.as_slice(), [1, 2]);
    }

    #[test]
    fn a_glob_matching_nothing_is_an_error() {
        assert!(matches!(eval("Nope*"), Err(NdxError::NoGroupMatches { .. })));
    }

    /// `atomid` works with nothing loaded — unlike every other predicate.
    #[test]
    fn atomid_needs_no_structure() {
        assert_eq!(eval("atomid 1-3").unwrap(), [1, 2, 3]);
        // And composes with the groups.
        assert_eq!(eval("atomid 1-5 & 1").unwrap(), [1, 2, 3]);
        assert_eq!(eval("1 & !atomid 1").unwrap(), [2, 3]);
    }

    #[test]
    fn unknown_group_is_reported() {
        assert!(matches!(
            eval("Nope"),
            Err(NdxError::UnknownGroup { .. })
        ));
    }

    /// Without -s / -p these still parse; the failure names the flag that would have worked.
    #[test]
    fn structure_predicates_say_what_is_missing() {
        match eval("element H") {
            Err(NdxError::NeedsSystem { feature, needs, hint, .. }) => {
                assert_eq!(feature, "element");
                assert_eq!(needs, "a structure file");
                assert!(hint.contains("-s conf.gro"));
            }
            other => panic!("expected NeedsSystem, got {other:?}"),
        }
        match eval("type OW") {
            Err(NdxError::NeedsSystem { feature, needs, hint, .. }) => {
                assert_eq!(feature, "type");
                assert_eq!(needs, "a topology");
                assert!(hint.contains("-p topol.top"));
            }
            other => panic!("expected NeedsSystem, got {other:?}"),
        }
    }

    #[test]
    fn bonded_says_it_needs_bonds() {
        match eval("element H & bonded 1") {
            Err(NdxError::NeedsSystem { feature, span, .. }) => {
                // The left operand is evaluated first, so it is the one that reports.
                assert_eq!(feature, "element");
                assert_eq!(span, 0..9, "the caret should sit under `element H`");
            }
            other => panic!("expected NeedsSystem, got {other:?}"),
        }
        match eval("bonded 1") {
            Err(NdxError::NeedsSystem { feature, needs, .. }) => {
                assert_eq!(feature, "bonded");
                assert_eq!(needs, "bond information");
            }
            other => panic!("expected NeedsSystem, got {other:?}"),
        }
    }

    #[test]
    fn within_says_it_needs_coordinates() {
        match eval("within 0.5 of 1") {
            Err(NdxError::NeedsSystem { feature, needs, .. }) => {
                assert_eq!(feature, "within");
                assert_eq!(needs, "atom coordinates");
            }
            other => panic!("expected NeedsSystem, got {other:?}"),
        }
    }

    #[test]
    fn rewrite_leaves_double_negation_alone() {
        // `!A & !B` has no universe-free form; it must stay an And.
        let e = rewrite(parse_expr("!1 & !2").unwrap(), true);
        assert!(matches!(e, Expr::And(..)));
    }

    #[test]
    fn rewrite_handles_either_side() {
        assert!(matches!(rewrite(parse_expr("1 & !2").unwrap(), true), Expr::Diff(..)));
        assert!(matches!(rewrite(parse_expr("!2 & 1").unwrap(), true), Expr::Diff(..)));
    }
}
