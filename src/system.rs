//! Loading the structure and topology that structure-dependent selections need.
//!
//! A [`SystemCtx`] is what turns `element H & bonded Protein` from an error into an answer. It is
//! optional: everything an `.ndx` can answer on its own still works with an empty one.

use std::path::PathBuf;

use crate::bonds::{self, BondGraph, BondSource};
use crate::error::{NdxError, Result};
use crate::model::IndexFile;
use crate::structure::Structure;
use crate::topology::Topology;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum BondMode {
    /// Real bonds from the topology; fall back to distances when there is no `.top`.
    #[default]
    Auto,
    /// Only the topology. Fails if none was given.
    Top,
    /// Always estimate from coordinates, even when a topology is available.
    Distance,
}

#[derive(Clone, Debug, Default)]
pub struct LoadOptions {
    pub structure: Option<PathBuf>,
    pub topology: Option<PathBuf>,
    /// `-D SYMBOL`, as in an `.mdp`'s `define = -DPOSRES`.
    pub defines: Vec<String>,
    /// `-I DIR` for `#include` resolution. `$GMXLIB` is searched too.
    pub include_dirs: Vec<PathBuf>,
    pub bonds: BondMode,
    pub pbc: bool,
}

impl LoadOptions {
    pub fn is_empty(&self) -> bool {
        self.structure.is_none() && self.topology.is_none()
    }
}

#[derive(Default)]
pub struct SystemCtx {
    pub structure: Option<Structure>,
    pub topology: Option<Topology>,
    pub bonds: Option<BondGraph>,
    pub pbc: bool,
    /// Things the user should know but that are not fatal.
    pub warnings: Vec<String>,
}

impl SystemCtx {
    pub fn is_empty(&self) -> bool {
        self.structure.is_none() && self.topology.is_none()
    }

    /// The number of atoms in the system, when we know it exactly.
    pub fn natoms(&self) -> Option<u32> {
        self.structure
            .as_ref()
            .map(Structure::natoms)
            .or_else(|| self.topology.as_ref().map(|t| t.natoms))
    }

    pub fn load(opts: &LoadOptions) -> Result<Self> {
        let mut cx = SystemCtx {
            pbc: opts.pbc,
            ..Default::default()
        };

        if let Some(p) = &opts.topology {
            let top = Topology::parse_path(p, &opts.defines, &opts.include_dirs)?;
            if !top.missing_includes.is_empty() {
                cx.warnings.push(format!(
                    "could not find {} #include(s): {}. \
                     That is fine if they only held force-field parameters; pass -I DIR (or set \
                     $GMXLIB) if a molecule definition is missing.",
                    top.missing_includes.len(),
                    top.missing_includes.join(", ")
                ));
            }
            if top.bonds.is_empty() {
                cx.warnings.push(format!(
                    "{} defines no bonds. If its water is rigid the bonds live in [ settles ], \
                     which is read; a coarse-grained model may genuinely have none.",
                    p.display()
                ));
            }
            cx.topology = Some(top);
        }

        if let Some(p) = &opts.structure {
            let mut s = Structure::from_gro_path(p)?;

            // A topology's masses identify elements far better than atom names do.
            if let Some(top) = &cx.topology {
                if top.natoms == s.natoms() {
                    s.refine_elements_from_masses(&top.masses);
                } else {
                    cx.warnings.push(format!(
                        "the structure has {} atoms but the topology expands to {}. \
                         They do not describe the same system; element and type lookups will be \
                         wrong. (Elements will be taken from atom names only.)",
                        s.natoms(),
                        top.natoms
                    ));
                }
            }

            if s.triclinic && opts.pbc {
                cx.warnings.push(
                    "the box is triclinic, but --pbc uses the rectangular minimum-image \
                     convention; distances across a boundary may be wrong"
                        .into(),
                );
            }

            cx.structure = Some(s);
        }

        cx.bonds = cx.build_bonds(opts)?;
        Ok(cx)
    }

    fn build_bonds(&mut self, opts: &LoadOptions) -> Result<Option<BondGraph>> {
        let from_top = |top: &Topology, natoms: u32| {
            BondGraph::from_pairs(natoms, &top.bonds, BondSource::Topology)
        };

        match opts.bonds {
            BondMode::Top => {
                let Some(top) = &self.topology else {
                    return Err(NdxError::Other(
                        "--bonds top needs a topology\nhint: pass -p topol.top".into(),
                    ));
                };
                let natoms = self.natoms().unwrap_or(top.natoms);
                Ok(Some(from_top(top, natoms)))
            }

            BondMode::Distance => {
                let Some(s) = &self.structure else {
                    return Err(NdxError::Other(
                        "--bonds distance needs coordinates\nhint: pass -s conf.gro".into(),
                    ));
                };
                Ok(Some(bonds::infer_from_distance(s)))
            }

            BondMode::Auto => match (&self.topology, &self.structure) {
                (Some(top), _) if !top.bonds.is_empty() => {
                    let natoms = self.natoms().unwrap_or(top.natoms);
                    Ok(Some(from_top(top, natoms)))
                }
                (_, Some(s)) => {
                    self.warnings.push(
                        "no topology bonds available, so bonds are being estimated from \
                         interatomic distances. Pass -p topol.top for the real connectivity."
                            .into(),
                    );
                    Ok(Some(bonds::infer_from_distance(s)))
                }
                _ => Ok(None),
            },
        }
    }

    /// Check the loaded system against the index file, so a mismatched pair is caught before it
    /// produces quietly wrong answers rather than after.
    pub fn check_against(&self, ndx: &IndexFile) -> Vec<String> {
        let mut out = Vec::new();
        let (Some(natoms), Some(max)) = (self.natoms(), ndx.max_atom()) else {
            return out;
        };
        if max > natoms {
            out.push(format!(
                "the index file references atom {max}, but the system has only {natoms} atoms. \
                 The .ndx and the structure/topology do not match."
            ));
        }
        out
    }
}

impl std::fmt::Debug for SystemCtx {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SystemCtx")
            .field("natoms", &self.natoms())
            .field("structure", &self.structure.is_some())
            .field("topology", &self.topology.is_some())
            .field("bonds", &self.bonds.as_ref().map(|b| (b.source, b.nbonds)))
            .finish()
    }
}

/// What a structure-dependent feature needs, for the error message when it is missing.
pub fn needs_structure(feature: &str) -> NdxError {
    NdxError::NeedsSystem {
        feature: feature.to_string(),
        needs: "a structure file",
        hint: "pass -s conf.gro",
        span: 0..0,
    }
}

pub fn needs_topology(feature: &str) -> NdxError {
    NdxError::NeedsSystem {
        feature: feature.to_string(),
        needs: "a topology",
        hint: "pass -p topol.top",
        span: 0..0,
    }
}

pub fn needs_bonds(feature: &str) -> NdxError {
    NdxError::NeedsSystem {
        feature: feature.to_string(),
        needs: "bond information",
        hint: "pass -p topol.top for real bonds, or -s conf.gro to estimate them from distances",
        span: 0..0,
    }
}

/// A summary line for `ndxed info` and the editor's banner.
pub fn describe(cx: &SystemCtx) -> Vec<String> {
    let mut out = Vec::new();
    if let Some(s) = &cx.structure {
        out.push(format!(
            "structure: {} atoms, box {:.3} x {:.3} x {:.3} nm{}",
            s.natoms(),
            s.box_diag[0],
            s.box_diag[1],
            s.box_diag[2],
            if s.triclinic { " (triclinic)" } else { "" }
        ));
    }
    if let Some(t) = &cx.topology {
        out.push(format!(
            "topology:  {} atoms, {} molecule type(s): {}",
            t.natoms,
            t.mol_names.len(),
            t.mol_names.join(", ")
        ));
    }
    if let Some(b) = &cx.bonds {
        out.push(format!(
            "bonds:     {} ({})",
            b.nbonds,
            match b.source {
                BondSource::Topology => "from the topology",
                BondSource::Distance => "estimated from distances",
            }
        ));
    }
    out
}
