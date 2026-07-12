use crate::error::{NdxError, Result};
use crate::model::{AtomId, Group, GroupRef, IndexFile};
use crate::system::{self, SystemCtx};

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum How {
    /// `K` near-equal parts. The first `len % K` parts get one extra atom.
    Parts(usize),
    /// Consecutive chunks of `K` atoms; the last one may be short.
    Size(usize),
    /// Cut at the given 1-based positions *within the group*.
    At(Vec<usize>),
    /// One group per residue — `make_ndx`'s `splitres`. Needs a structure.
    Residue,
    /// One group per molecule instance — the GROMACS answer to `splitch`. Needs a topology.
    Molecule,
}

pub struct SplitResult {
    pub names: Vec<String>,
    pub sizes: Vec<usize>,
}

/// Split a group into several groups. Chunks follow the group's own atom order, so splitting an
/// unsorted group does not silently reorder it.
pub fn split(
    ndx: &mut IndexFile,
    r: &GroupRef,
    how: &How,
    prefix: Option<&str>,
    replace: bool,
    system: &SystemCtx,
) -> Result<SplitResult> {
    let id = ndx.resolve(r)?;
    let src = ndx.get(id)?;
    let prefix = prefix.unwrap_or(&src.name).to_string();

    // Each chunk carries a label; the positional splits have none and get numbered.
    let chunks: Vec<(Option<String>, Vec<AtomId>)> = match how {
        How::Residue | How::Molecule => by_key(&src.atoms, how, system)?,
        _ => chunk(&src.atoms, how)?
            .into_iter()
            .map(|atoms| (None, atoms))
            .collect(),
    };

    let mut names = Vec::with_capacity(chunks.len());
    let mut sizes = Vec::with_capacity(chunks.len());

    // Build the groups before touching `ndx`, so a failure leaves it untouched.
    let mut new_groups = Vec::with_capacity(chunks.len());
    for (i, (label, atoms)) in chunks.into_iter().enumerate() {
        let name = match label {
            Some(l) => format!("{prefix}_{l}"),
            None => format!("{prefix}_{}", i + 1),
        };
        sizes.push(atoms.len());
        names.push(name.clone());
        new_groups.push(Group::new(name, atoms));
    }

    if replace {
        ndx.remove_many(&[id]);
    }
    for g in new_groups {
        ndx.push(g);
    }

    Ok(SplitResult { names, sizes })
}

/// One group per residue or per molecule instance.
///
/// The parts come out in first-appearance order and each keeps the source group's atom order, so a
/// group whose atoms are interleaved across residues still splits cleanly.
fn by_key(
    atoms: &[AtomId],
    how: &How,
    system: &SystemCtx,
) -> Result<Vec<(Option<String>, Vec<AtomId>)>> {
    // (sort key, label) per atom.
    let keyed: Vec<(i64, String)> = match how {
        How::Residue => {
            let s = system
                .structure
                .as_ref()
                .ok_or_else(|| system::needs_structure("split --by residue"))?;
            atoms
                .iter()
                .map(|&a| {
                    let resid = s.resid(a).unwrap_or(0);
                    let resname = s.resname(a).unwrap_or("UNK");
                    (resid as i64, format!("{resname}{resid}"))
                })
                .collect()
        }
        How::Molecule => {
            let t = system
                .topology
                .as_ref()
                .ok_or_else(|| system::needs_topology("split --by molecule"))?;
            atoms
                .iter()
                .map(|&a| {
                    let inst = t.instance(a).unwrap_or(0);
                    let name = t.molname(a).unwrap_or("UNK");
                    // Distinguish the copies of a species, and different species from each other.
                    (
                        (t.mol_names.iter().position(|n| n == name).unwrap_or(0) as i64) << 32
                            | inst as i64,
                        format!("{name}{}", inst + 1),
                    )
                })
                .collect()
        }
        _ => unreachable!("by_key is only for Residue/Molecule"),
    };

    let mut order: Vec<i64> = Vec::new();
    let mut buckets: std::collections::HashMap<i64, (String, Vec<AtomId>)> =
        std::collections::HashMap::new();
    for (&a, (key, label)) in atoms.iter().zip(keyed) {
        let e = buckets.entry(key).or_insert_with(|| {
            order.push(key);
            (label, Vec::new())
        });
        e.1.push(a);
    }

    Ok(order
        .into_iter()
        .map(|k| {
            let (label, atoms) = buckets.remove(&k).expect("key came from the map");
            (Some(label), atoms)
        })
        .collect())
}

fn chunk(atoms: &[AtomId], how: &How) -> Result<Vec<Vec<AtomId>>> {
    let len = atoms.len();
    match how {
        How::Parts(k) => {
            let k = *k;
            if k == 0 {
                return Err(NdxError::Other("--parts must be at least 1".into()));
            }
            if k > len {
                return Err(NdxError::SplitTooFine { parts: k, len });
            }
            // len = q*k + r: the first r parts take q+1 atoms, the rest take q.
            let (q, r) = (len / k, len % k);
            let mut out = Vec::with_capacity(k);
            let mut i = 0;
            for p in 0..k {
                let take = q + usize::from(p < r);
                out.push(atoms[i..i + take].to_vec());
                i += take;
            }
            Ok(out)
        }

        How::Size(k) => {
            let k = *k;
            if k == 0 {
                return Err(NdxError::Other("--size must be at least 1".into()));
            }
            Ok(atoms.chunks(k).map(<[AtomId]>::to_vec).collect())
        }

        How::Residue | How::Molecule => unreachable!("handled by by_key"),

        How::At(cuts) => {
            let mut cuts: Vec<usize> = cuts.clone();
            cuts.sort_unstable();
            cuts.dedup();
            for &c in &cuts {
                if c == 0 || c > len {
                    return Err(NdxError::SplitBoundaryOutOfRange { at: c, len });
                }
            }
            let mut out = Vec::with_capacity(cuts.len() + 1);
            let mut start = 0usize;
            for &c in &cuts {
                // `--at 100` means "cut after the 100th atom", so the first part is 1..=100.
                if c > start {
                    out.push(atoms[start..c].to_vec());
                    start = c;
                }
            }
            if start < len {
                out.push(atoms[start..].to_vec());
            }
            Ok(out)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ndx() -> IndexFile {
        IndexFile {
            groups: vec![Group::new("P", vec![10, 20, 30, 40, 50, 60, 70])],
        }
    }

    fn sizes(how: How) -> Vec<usize> {
        let mut f = ndx();
        split(&mut f, &GroupRef::Id(0), &how, None, false, &SystemCtx::default())
            .unwrap()
            .sizes
    }

    #[test]
    fn parts_distributes_the_remainder_to_the_front() {
        // 7 atoms into 3 parts -> 3, 2, 2
        assert_eq!(sizes(How::Parts(3)), [3, 2, 2]);
    }

    #[test]
    fn parts_that_divide_evenly() {
        let mut f = IndexFile {
            groups: vec![Group::new("P", (1..=6).collect())],
        };
        let r = split(&mut f, &GroupRef::Id(0), &How::Parts(3), None, false, &SystemCtx::default()).unwrap();
        assert_eq!(r.sizes, [2, 2, 2]);
    }

    #[test]
    fn size_leaves_a_short_tail() {
        assert_eq!(sizes(How::Size(3)), [3, 3, 1]);
    }

    #[test]
    fn at_cuts_after_the_given_positions() {
        assert_eq!(sizes(How::At(vec![2, 5])), [2, 3, 2]);
    }

    #[test]
    fn chunks_follow_file_order() {
        let mut f = IndexFile {
            groups: vec![Group::new("P", vec![5, 1, 4, 2, 3])],
        };
        split(&mut f, &GroupRef::Id(0), &How::Size(2), None, false, &SystemCtx::default()).unwrap();
        assert_eq!(f.groups[1].atoms, [5, 1]);
        assert_eq!(f.groups[2].atoms, [4, 2]);
        assert_eq!(f.groups[3].atoms, [3]);
    }

    #[test]
    fn names_are_prefixed_and_numbered() {
        let mut f = ndx();
        let r = split(&mut f, &GroupRef::Id(0), &How::Parts(2), None, false, &SystemCtx::default()).unwrap();
        assert_eq!(r.names, ["P_1", "P_2"]);
    }

    #[test]
    fn custom_prefix() {
        let mut f = ndx();
        let r = split(&mut f, &GroupRef::Id(0), &How::Parts(2), Some("chunk"), false, &SystemCtx::default()).unwrap();
        assert_eq!(r.names, ["chunk_1", "chunk_2"]);
    }

    #[test]
    fn replace_drops_the_source() {
        let mut f = ndx();
        split(&mut f, &GroupRef::Id(0), &How::Parts(2), None, true, &SystemCtx::default()).unwrap();
        let names: Vec<&str> = f.groups.iter().map(|g| g.name.as_str()).collect();
        assert_eq!(names, ["P_1", "P_2"]);
    }

    #[test]
    fn without_replace_the_source_stays() {
        let mut f = ndx();
        split(&mut f, &GroupRef::Id(0), &How::Parts(2), None, false, &SystemCtx::default()).unwrap();
        assert_eq!(f.len(), 3);
        assert_eq!(f.groups[0].name, "P");
    }

    #[test]
    fn too_many_parts_is_an_error() {
        let mut f = ndx();
        assert!(matches!(
            split(&mut f, &GroupRef::Id(0), &How::Parts(99), None, false, &SystemCtx::default()),
            Err(NdxError::SplitTooFine { parts: 99, len: 7 })
        ));
        assert_eq!(f.len(), 1, "the file must be untouched on failure");
    }

    #[test]
    fn out_of_range_boundary_is_an_error() {
        let mut f = ndx();
        assert!(matches!(
            split(&mut f, &GroupRef::Id(0), &How::At(vec![99]), None, false, &SystemCtx::default()),
            Err(NdxError::SplitBoundaryOutOfRange { at: 99, len: 7 })
        ));
        assert!(matches!(
            split(&mut f, &GroupRef::Id(0), &How::At(vec![0]), None, false, &SystemCtx::default()),
            Err(NdxError::SplitBoundaryOutOfRange { at: 0, len: 7 })
        ));
    }

    #[test]
    fn a_cut_at_the_very_end_yields_one_part() {
        assert_eq!(sizes(How::At(vec![7])), [7]);
    }

    #[test]
    fn every_atom_survives_a_split() {
        for how in [How::Parts(3), How::Size(3), How::At(vec![2, 5])] {
            let mut f = ndx();
            split(&mut f, &GroupRef::Id(0), &how, None, true, &SystemCtx::default()).unwrap();
            let all: Vec<AtomId> = f.groups.iter().flat_map(|g| g.atoms.clone()).collect();
            assert_eq!(all, [10, 20, 30, 40, 50, 60, 70], "for {how:?}");
        }
    }
}
