//! What `!` complements against.
//!
//! Without a structure file the set of all atoms is not knowable from the index file alone, so we
//! derive it and say so. Note that the most common `!` idiom, `A & !B`, is rewritten to a plain
//! difference in `expr::eval` and never gets here — see the `rewrite` function there.

use crate::atomset::AtomSet;
use crate::error::{NdxError, Result};
use crate::model::{GroupRef, IndexFile};
use crate::system::SystemCtx;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum UniverseSpec {
    /// Derive it: structure file, else `System`, else the union of all groups.
    #[default]
    Auto,
    Natoms(u32),
    Group(GroupRef),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum UniverseOrigin {
    /// From a loaded structure file. Exact.
    Structure(u32),
    /// From `--natoms`.
    Explicit(u32),
    /// From `--universe GROUP`, or the file's own `System` group.
    NamedGroup { id: usize, name: String },
    /// Nothing better was available.
    UnionOfAll,
}

#[derive(Clone, Debug)]
pub struct Universe {
    pub set: AtomSet,
    pub origin: UniverseOrigin,
    /// Shown once on stderr when the universe had to be guessed.
    pub warning: Option<String>,
}

pub fn resolve(ndx: &IndexFile, spec: &UniverseSpec, system: &SystemCtx) -> Result<Universe> {
    let union_all = ndx.union_all();

    // A loaded structure or topology settles it exactly, so `!` needs no guessing at all.
    if let Some(natoms) = system.natoms() {
        let warning = union_all.max().filter(|m| *m > natoms).map(|m| {
            format!(
                "the index file references atom {m}, but the loaded system has only {natoms} \
                 atoms; '!' will not see the extra atoms"
            )
        });
        return Ok(Universe {
            set: AtomSet::range_inclusive(1, natoms),
            origin: UniverseOrigin::Structure(natoms),
            warning,
        });
    }

    match spec {
        UniverseSpec::Natoms(n) => {
            if let Some(max) = union_all.max()
                && max > *n
            {
                return Err(NdxError::AtomExceedsNatoms {
                    idx: max,
                    natoms: *n,
                });
            }
            Ok(Universe {
                set: AtomSet::range_inclusive(1, *n),
                origin: UniverseOrigin::Explicit(*n),
                warning: None,
            })
        }

        UniverseSpec::Group(r) => {
            let id = ndx.resolve(r)?;
            let g = ndx.get(id)?;
            let set = g.to_set();
            let warning = (!union_all.is_subset_of(&set)).then(|| {
                format!(
                    "the universe group {:?} does not contain every atom in the file \
                     ({} atom(s) fall outside it); those atoms cannot appear in any complement",
                    g.name,
                    union_all.difference(&set).len()
                )
            });
            Ok(Universe {
                set,
                origin: UniverseOrigin::NamedGroup {
                    id,
                    name: g.name.clone(),
                },
                warning,
            })
        }

        UniverseSpec::Auto => {
            if ndx.is_empty() {
                return Err(NdxError::NoUniverse);
            }
            if let Some(id) = find_system(ndx) {
                let g = ndx.get(id)?;
                let sys = g.to_set();
                // A `System` group that doesn't cover the file is corrupt. Widening is the safe
                // reading: a too-small universe silently drops atoms from every complement.
                if union_all.is_subset_of(&sys) {
                    return Ok(Universe {
                        set: sys,
                        origin: UniverseOrigin::NamedGroup {
                            id,
                            name: g.name.clone(),
                        },
                        warning: None,
                    });
                }
                let missing = union_all.difference(&sys).len();
                return Ok(Universe {
                    set: sys.union(&union_all),
                    origin: UniverseOrigin::NamedGroup {
                        id,
                        name: g.name.clone(),
                    },
                    warning: Some(format!(
                        "the {:?} group is missing {missing} atom(s) that other groups reference; \
                         using its union with all groups as the universe for '!'",
                        g.name
                    )),
                });
            }

            let n = union_all.len();
            let max = union_all.max().unwrap_or(0);
            Ok(Universe {
                set: union_all,
                origin: UniverseOrigin::UnionOfAll,
                warning: Some(format!(
                    "'!' is complementing against the union of all groups ({n} atoms, highest \
                     index {max}) because the file has no 'System' group. Atoms that exist in the \
                     simulation but appear in no group will be missing from the result.\n\
                     hint: pass --natoms N (or --universe GROUP) to be exact"
                )),
            })
        }
    }
}

/// The one group named `System` (case-insensitive fallback). Ambiguous duplicates are ignored —
/// we would rather fall through to the union than silently pick one of two `System` groups.
fn find_system(ndx: &IndexFile) -> Option<usize> {
    let exact = ndx.find_by_name("System");
    if exact.len() == 1 {
        return Some(exact[0]);
    }
    if exact.is_empty() {
        let ci: Vec<usize> = ndx
            .groups
            .iter()
            .enumerate()
            .filter(|(_, g)| g.name.eq_ignore_ascii_case("system"))
            .map(|(i, _)| i)
            .collect();
        if ci.len() == 1 {
            return Some(ci[0]);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Group;

    fn ndx(groups: Vec<Group>) -> IndexFile {
        IndexFile { groups }
    }

    #[test]
    fn system_group_wins() {
        let f = ndx(vec![
            Group::new("System", (1..=10).collect()),
            Group::new("Protein", vec![1, 2]),
        ]);
        let u = resolve(&f, &UniverseSpec::Auto, &SystemCtx::default()).unwrap();
        assert_eq!(u.set.len(), 10);
        assert!(u.warning.is_none());
        assert!(matches!(u.origin, UniverseOrigin::NamedGroup { id: 0, .. }));
    }

    #[test]
    fn system_matched_case_insensitively() {
        let f = ndx(vec![Group::new("SYSTEM", (1..=4).collect())]);
        let u = resolve(&f, &UniverseSpec::Auto, &SystemCtx::default()).unwrap();
        assert!(matches!(u.origin, UniverseOrigin::NamedGroup { .. }));
    }

    #[test]
    fn without_system_falls_back_to_union_and_warns() {
        let f = ndx(vec![
            Group::new("Protein", vec![1, 2]),
            Group::new("SOL", vec![3, 4]),
        ]);
        let u = resolve(&f, &UniverseSpec::Auto, &SystemCtx::default()).unwrap();
        assert_eq!(u.set.as_slice(), [1, 2, 3, 4]);
        assert_eq!(u.origin, UniverseOrigin::UnionOfAll);
        assert!(u.warning.is_some());
    }

    #[test]
    fn undersized_system_is_widened_with_a_warning() {
        let f = ndx(vec![
            Group::new("System", vec![1, 2]),
            Group::new("Extra", vec![3]),
        ]);
        let u = resolve(&f, &UniverseSpec::Auto, &SystemCtx::default()).unwrap();
        assert_eq!(u.set.as_slice(), [1, 2, 3]);
        assert!(u.warning.is_some());
    }

    #[test]
    fn duplicate_system_groups_fall_through() {
        let f = ndx(vec![
            Group::new("System", vec![1]),
            Group::new("System", vec![2]),
        ]);
        let u = resolve(&f, &UniverseSpec::Auto, &SystemCtx::default()).unwrap();
        assert_eq!(u.origin, UniverseOrigin::UnionOfAll);
    }

    #[test]
    fn natoms_is_exact_and_silent() {
        let f = ndx(vec![Group::new("Protein", vec![1, 2])]);
        let u = resolve(&f, &UniverseSpec::Natoms(100), &SystemCtx::default()).unwrap();
        assert_eq!(u.set.len(), 100);
        assert!(u.warning.is_none());
    }

    #[test]
    fn natoms_smaller_than_the_file_is_an_error() {
        let f = ndx(vec![Group::new("Protein", vec![1, 50])]);
        assert!(matches!(
            resolve(&f, &UniverseSpec::Natoms(10), &SystemCtx::default()),
            Err(NdxError::AtomExceedsNatoms { idx: 50, natoms: 10 })
        ));
    }

    #[test]
    fn explicit_universe_group() {
        let f = ndx(vec![
            Group::new("All", (1..=5).collect()),
            Group::new("A", vec![1]),
        ]);
        let spec = UniverseSpec::Group(GroupRef::name("All"));
        let u = resolve(&f, &spec, &SystemCtx::default()).unwrap();
        assert_eq!(u.set.len(), 5);
        assert!(u.warning.is_none());
    }

    #[test]
    fn explicit_universe_that_misses_atoms_warns() {
        let f = ndx(vec![
            Group::new("Small", vec![1]),
            Group::new("A", vec![1, 2, 3]),
        ]);
        let spec = UniverseSpec::Group(GroupRef::name("Small"));
        let u = resolve(&f, &spec, &SystemCtx::default()).unwrap();
        assert!(u.warning.is_some());
    }

    #[test]
    fn no_groups_no_universe() {
        assert!(matches!(
            resolve(&IndexFile::new(), &UniverseSpec::Auto, &SystemCtx::default()),
            Err(NdxError::NoUniverse)
        ));
    }
}
