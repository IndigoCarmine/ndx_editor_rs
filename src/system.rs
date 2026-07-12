//! Structure and topology: the extension point for `.gro` / `.top` support.
//!
//! Nothing here is implemented yet — [`SystemCtx`] is always empty in the current build, and the
//! expression evaluator reports a specific "you need `-s conf.gro`" error when an expression
//! reaches for something a `SystemCtx` would provide.
//!
//! The point of defining the traits now is that the expression grammar, the auto-naming, the CLI
//! shape and the tests around them are all already written against them. Adding `.gro` and `.top`
//! is then purely additive: implement a reader, fill in the matching arm in `expr::eval`, add the
//! flag. No existing code changes.
//!
//! Planned implementations:
//!
//! - `GroStructure` / `PdbStructure`  -> [`Structure`]
//! - `TopTopology` (`#include` expansion, `[moleculetype]`, `[molecules]` expansion) -> [`Topology`]
//! - `TopBonds` (from `[bonds]`, the real bond graph)                                -> [`BondGraph`]
//! - `DistanceBonds` (coordinates + covalent radii, cell list; the fallback when no `.top`)

use crate::model::AtomId;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Element {
    H,
    C,
    N,
    O,
    P,
    S,
    Other,
}

/// Per-atom identity and coordinates, from a `.gro` or `.pdb`.
pub trait Structure {
    fn natoms(&self) -> u32;
    fn name(&self, a: AtomId) -> &str;
    fn resname(&self, a: AtomId) -> &str;
    fn resid(&self, a: AtomId) -> i32;
    fn chain(&self, a: AtomId) -> Option<char>;
    /// `.gro` has no element column, so this is inferred from the atom name.
    fn element(&self, a: AtomId) -> Option<Element>;
    /// Position in nm.
    fn pos(&self, a: AtomId) -> [f32; 3];
}

/// Per-atom force-field properties, from a `.top` / `.itp`.
pub trait Topology {
    fn atomtype(&self, a: AtomId) -> &str;
    fn charge(&self, a: AtomId) -> f32;
    fn mass(&self, a: AtomId) -> f32;
}

/// What "adjacent" means. Two implementations are planned: the real bond list out of a `.top`,
/// and a distance-based guess from `.gro` coordinates when no topology is available.
pub trait BondGraph {
    fn neighbors(&self, a: AtomId) -> &[AtomId];
}

/// Everything we know about the system beyond the index file itself.
///
/// Always empty today. `EvalCtx` and `Session` already take one, so wiring up `-s` / `-p` later
/// touches only the CLI and the readers.
#[derive(Default)]
pub struct SystemCtx {
    pub structure: Option<Box<dyn Structure>>,
    pub topology: Option<Box<dyn Topology>>,
    pub bonds: Option<Box<dyn BondGraph>>,
}

impl SystemCtx {
    pub fn is_empty(&self) -> bool {
        self.structure.is_none() && self.topology.is_none() && self.bonds.is_none()
    }
}

impl std::fmt::Debug for SystemCtx {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SystemCtx")
            .field("structure", &self.structure.is_some())
            .field("topology", &self.topology.is_some())
            .field("bonds", &self.bonds.is_some())
            .finish()
    }
}
