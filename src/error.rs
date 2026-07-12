use std::ops::Range;
use std::path::PathBuf;

use crate::model::AtomId;

/// Process exit codes. `diff` uses 1 for "files differ", like diff(1).
pub mod exit {
    pub const OK: i32 = 0;
    pub const DIFFERENT: i32 = 1;
    pub const USAGE: i32 = 2;
    pub const DOMAIN: i32 = 3;
    pub const IO: i32 = 4;
}

pub type Result<T> = std::result::Result<T, NdxError>;

#[derive(Debug, thiserror::Error)]
pub enum NdxError {
    #[error("{path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error(transparent)]
    Parse(#[from] ParseError),

    #[error(transparent)]
    Expr(#[from] ExprError),

    #[error("no group named {name:?}{}", suggestion_suffix(.suggestion))]
    UnknownGroup {
        name: String,
        suggestion: Option<String>,
    },

    #[error("group id {id} is out of range: the file has {count} group(s)")]
    GroupIdOutOfRange { id: usize, count: usize },

    #[error(
        "group name {name:?} is ambiguous: it matches ids {ids:?}\n\
         hint: use a numeric id, or {name}#1 / {name}#2 to pick an occurrence"
    )]
    AmbiguousName { name: String, ids: Vec<usize> },

    #[error("group {name:?} has only {have} occurrence(s), but #{want} was requested")]
    NoSuchOccurrence {
        name: String,
        have: usize,
        want: usize,
    },

    #[error("no group name matches the pattern {pattern:?}")]
    NoGroupMatches { pattern: String },

    #[error(
        "{pattern:?} matches {} groups ({}), but this takes exactly one\n\
         hint: name one of them, or use its id",
        .names.len(),
        .names.join(", ")
    )]
    PatternNotUnique {
        pattern: String,
        names: Vec<String>,
    },

    #[error(
        "cannot evaluate '!' (complement): the index file has no groups, so there is no \
         universe of atoms to complement against\n\
         hint: pass --natoms N"
    )]
    NoUniverse,

    #[error("atom index {idx} exceeds --natoms {natoms}")]
    AtomExceedsNatoms { idx: AtomId, natoms: u32 },

    #[error("the expression produced an empty group")]
    EmptyResult,

    #[error("--parts {parts} exceeds the group's {len} atom(s)")]
    SplitTooFine { parts: usize, len: usize },

    #[error("split boundary {at} is outside the group's 1..={len} atoms")]
    SplitBoundaryOutOfRange { at: usize, len: usize },

    #[error("group {group:?} has only {len} atom(s), but {n} were requested")]
    NotEnoughAtoms {
        group: String,
        len: usize,
        n: usize,
    },

    /// A structure-dependent feature was used without a structure file.
    /// The gro/top readers are not implemented yet (M10/M11).
    #[error(
        "`{feature}` needs {needs}, which has not been loaded\n\
         hint: {hint}"
    )]
    NeedsSystem {
        feature: String,
        needs: &'static str,
        hint: &'static str,
        span: Range<usize>,
    },

    #[error("refusing to overwrite the input file {path}\n\
             hint: write somewhere else, or pass --force-overwrite")]
    WouldOverwriteInput { path: PathBuf },

    #[error("no output file set\nhint: use `q FILE` / `w FILE`, or start with -o FILE")]
    NoOutputPath,

    #[error("{0}")]
    Other(String),
}

fn suggestion_suffix(s: &Option<String>) -> String {
    match s {
        Some(s) => format!("\nhint: did you mean {s:?}?"),
        None => String::new(),
    }
}

impl NdxError {
    pub fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        NdxError::Io {
            path: path.into(),
            source,
        }
    }

    pub fn exit_code(&self) -> i32 {
        match self {
            NdxError::Io { .. } => exit::IO,
            _ => exit::DOMAIN,
        }
    }

    /// The source span within the user's expression, if this error came from one.
    pub fn span(&self) -> Option<Range<usize>> {
        match self {
            NdxError::Expr(e) => e.span(),
            NdxError::NeedsSystem { span, .. } => Some(span.clone()),
            _ => None,
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum ParseError {
    #[error("{origin}:{line}: unterminated group header (expected ']')")]
    UnterminatedHeader { origin: String, line: usize },

    #[error("{origin}:{line}: empty group name")]
    EmptyName { origin: String, line: usize },

    #[error("{origin}:{line}:{col}: {token:?} is not an atom index")]
    BadToken {
        origin: String,
        line: usize,
        col: usize,
        token: String,
    },

    #[error("{origin}:{line}:{col}: atom index 0 is invalid (.ndx indices are 1-based)")]
    ZeroIndex {
        origin: String,
        line: usize,
        col: usize,
    },

    #[error("{origin}:{line}:{col}: atom indices appear before any '[ group ]' header")]
    AtomsBeforeGroup {
        origin: String,
        line: usize,
        col: usize,
    },

    #[error("{origin}:{line}:{col}: atom index {token} is too large")]
    IndexTooLarge {
        origin: String,
        line: usize,
        col: usize,
        token: String,
    },
}

#[derive(Debug, thiserror::Error)]
pub enum ExprError {
    #[error("empty expression")]
    Empty,

    #[error("expected a group reference")]
    ExpectedGroupRef { span: Range<usize> },

    #[error("unbalanced parenthesis: expected ')'")]
    UnbalancedParen { span: Range<usize> },

    #[error("unterminated quoted name")]
    UnterminatedQuote { span: Range<usize> },

    #[error("unexpected trailing input")]
    TrailingInput { span: Range<usize> },

    #[error("'{op}' is missing an operand")]
    MissingOperand { op: String, span: Range<usize> },

    #[error("`{keyword}` needs an argument")]
    MissingArgument {
        keyword: &'static str,
        span: Range<usize>,
    },

    #[error("expected `of` after `within <radius>`")]
    ExpectedOf { span: Range<usize> },

    #[error("{msg}")]
    BadArgument { msg: String, span: Range<usize> },
}

impl ExprError {
    pub fn span(&self) -> Option<Range<usize>> {
        match self {
            ExprError::Empty => None,
            ExprError::ExpectedGroupRef { span }
            | ExprError::UnbalancedParen { span }
            | ExprError::UnterminatedQuote { span }
            | ExprError::TrailingInput { span }
            | ExprError::MissingOperand { span, .. }
            | ExprError::MissingArgument { span, .. }
            | ExprError::ExpectedOf { span }
            | ExprError::BadArgument { span, .. } => Some(span.clone()),
        }
    }
}

/// Render a caret diagnostic under the offending part of an expression:
///
/// ```text
///   element H & bonded Protein
///               ^^^^^^^^^^^^^^
/// ```
pub fn render_caret(src: &str, span: &Range<usize>) -> String {
    let start = span.start.min(src.len());
    let end = span.end.clamp(start, src.len());
    // Column in characters, so multi-byte group names line up.
    let pad = src[..start].chars().count();
    let width = src[start..end].chars().count().max(1);
    format!("  {src}\n  {}{}", " ".repeat(pad), "^".repeat(width))
}
