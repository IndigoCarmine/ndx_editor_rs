pub mod diff;
pub mod merge;
pub mod select;
pub mod split;

use crate::atomset::AtomSet;
use crate::error::{NdxError, Result};
use crate::model::{Group, GroupRef, IndexFile};

/// Rename a group in place.
pub fn rename(ndx: &mut IndexFile, r: &GroupRef, new_name: &str) -> Result<String> {
    let id = ndx.resolve(r)?;
    let g = ndx.get_mut(id)?;
    let old = std::mem::replace(&mut g.name, new_name.to_string());
    Ok(old)
}

/// Delete groups. Remaining ids shift down.
pub fn delete(ndx: &mut IndexFile, refs: &[GroupRef]) -> Result<Vec<String>> {
    let ids = ndx.resolve_many(refs)?;
    let names = ids.iter().map(|&i| ndx.groups[i].name.clone()).collect();
    ndx.remove_many(&ids);
    Ok(names)
}

/// Keep only these groups, in the order given.
pub fn keep(ndx: &mut IndexFile, refs: &[GroupRef]) -> Result<()> {
    let ids = ndx.resolve_many(refs)?;
    ndx.retain_ids(&ids);
    Ok(())
}

/// A new group from literal atom numbers: `1-10,15,20-30`.
pub fn atoms(ndx: &mut IndexFile, ranges: &str, name: Option<&str>) -> Result<usize> {
    let set = crate::atomset::parse_ranges(ranges)?;
    let name = match name {
        Some(n) => n.to_string(),
        None => ndx.unique_name(&format!("atoms_{ranges}")),
    };
    Ok(ndx.push(Group::new(name, set.into_vec())))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum End {
    Head,
    Tail,
}

/// The first (or last) `n` atoms of a group, **in file order** — "first 100" means the first 100
/// as written, not the 100 lowest indices.
pub fn head(
    ndx: &mut IndexFile,
    r: &GroupRef,
    n: usize,
    end: End,
    name: Option<&str>,
    strict: bool,
) -> Result<usize> {
    let id = ndx.resolve(r)?;
    let src = ndx.get(id)?;
    if strict && n > src.len() {
        return Err(NdxError::NotEnoughAtoms {
            group: src.name.clone(),
            len: src.len(),
            n,
        });
    }
    let take = n.min(src.len());
    let atoms = match end {
        End::Head => src.atoms[..take].to_vec(),
        End::Tail => src.atoms[src.len() - take..].to_vec(),
    };
    let default = match end {
        End::Head => format!("{}_head{n}", src.name),
        End::Tail => format!("{}_tail{n}", src.name),
    };
    let name = name.map(str::to_owned).unwrap_or_else(|| ndx.unique_name(&default));
    Ok(ndx.push(Group::new(name, atoms)))
}

/// Canonicalize groups: sort and/or drop duplicate atoms.
pub fn fmt(ndx: &mut IndexFile, sort: bool, dedup: bool) {
    for g in &mut ndx.groups {
        if sort && dedup {
            g.atoms = AtomSet::from_unsorted(std::mem::take(&mut g.atoms)).into_vec();
        } else if sort {
            g.atoms.sort_unstable();
        } else if dedup {
            let mut seen = std::collections::HashSet::new();
            g.atoms.retain(|a| seen.insert(*a));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ndx() -> IndexFile {
        IndexFile {
            groups: vec![
                Group::new("System", (1..=10).collect()),
                // Deliberately unsorted: file order must survive head/tail/split.
                Group::new("Protein", vec![5, 1, 4, 2, 3]),
            ],
        }
    }

    #[test]
    fn rename_returns_the_old_name() {
        let mut f = ndx();
        assert_eq!(rename(&mut f, &GroupRef::Id(1), "Prot").unwrap(), "Protein");
        assert_eq!(f.groups[1].name, "Prot");
    }

    #[test]
    fn delete_shifts_ids() {
        let mut f = ndx();
        delete(&mut f, &[GroupRef::Id(0)]).unwrap();
        assert_eq!(f.len(), 1);
        assert_eq!(f.groups[0].name, "Protein");
    }

    #[test]
    fn keep_reorders() {
        let mut f = ndx();
        keep(&mut f, &[GroupRef::Id(1), GroupRef::Id(0)]).unwrap();
        let names: Vec<&str> = f.groups.iter().map(|g| g.name.as_str()).collect();
        assert_eq!(names, ["Protein", "System"]);
    }

    #[test]
    fn head_uses_file_order_not_sorted_order() {
        let mut f = ndx();
        let id = head(&mut f, &GroupRef::Id(1), 3, End::Head, None, false).unwrap();
        assert_eq!(f.groups[id].atoms, [5, 1, 4]);
        assert_eq!(f.groups[id].name, "Protein_head3");
    }

    #[test]
    fn tail_uses_file_order() {
        let mut f = ndx();
        let id = head(&mut f, &GroupRef::Id(1), 2, End::Tail, None, false).unwrap();
        assert_eq!(f.groups[id].atoms, [2, 3]);
    }

    #[test]
    fn head_clamps_by_default() {
        let mut f = ndx();
        let id = head(&mut f, &GroupRef::Id(1), 99, End::Head, None, false).unwrap();
        assert_eq!(f.groups[id].atoms.len(), 5);
    }

    #[test]
    fn head_strict_errors() {
        let mut f = ndx();
        assert!(matches!(
            head(&mut f, &GroupRef::Id(1), 99, End::Head, None, true),
            Err(NdxError::NotEnoughAtoms { len: 5, n: 99, .. })
        ));
    }

    #[test]
    fn head_of_zero_is_empty() {
        let mut f = ndx();
        let id = head(&mut f, &GroupRef::Id(1), 0, End::Head, None, false).unwrap();
        assert!(f.groups[id].atoms.is_empty());
    }

    #[test]
    fn atoms_from_ranges() {
        let mut f = ndx();
        let id = atoms(&mut f, "1-3,7", Some("Mine")).unwrap();
        assert_eq!(f.groups[id].atoms, [1, 2, 3, 7]);
        assert_eq!(f.groups[id].name, "Mine");
    }

    #[test]
    fn fmt_sorts_and_dedups() {
        let mut f = IndexFile {
            groups: vec![Group::new("A", vec![3, 1, 1, 2])],
        };
        fmt(&mut f, true, true);
        assert_eq!(f.groups[0].atoms, [1, 2, 3]);
    }

    #[test]
    fn fmt_dedup_alone_keeps_order() {
        let mut f = IndexFile {
            groups: vec![Group::new("A", vec![3, 1, 3, 2])],
        };
        fmt(&mut f, false, true);
        assert_eq!(f.groups[0].atoms, [3, 1, 2]);
    }
}
