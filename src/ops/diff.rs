use crate::atomset::AtomSet;
use crate::model::IndexFile;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum Match {
    /// Pair groups up by name.
    #[default]
    Name,
    /// Pair groups up by position.
    Index,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Change {
    Added,
    Removed,
    /// Same name, same atoms, but written in a different order (or with different duplicates).
    Reordered,
    Changed {
        added: AtomSet,
        removed: AtomSet,
    },
    Same,
}

#[derive(Clone, Debug)]
pub struct Entry {
    pub name: String,
    pub change: Change,
    pub len_a: Option<usize>,
    pub len_b: Option<usize>,
}

pub struct Diff {
    pub entries: Vec<Entry>,
}

impl Diff {
    pub fn is_empty(&self) -> bool {
        self.entries.iter().all(|e| e.change == Change::Same)
    }
}

pub fn diff(a: &IndexFile, b: &IndexFile, by: Match) -> Diff {
    match by {
        Match::Index => diff_by_index(a, b),
        Match::Name => diff_by_name(a, b),
    }
}

fn diff_by_index(a: &IndexFile, b: &IndexFile) -> Diff {
    let mut entries = Vec::new();
    for i in 0..a.len().max(b.len()) {
        match (a.groups.get(i), b.groups.get(i)) {
            (Some(ga), Some(gb)) => {
                // A rename at the same position reads as a change, which is what --by index is for.
                let name = if ga.name == gb.name {
                    ga.name.clone()
                } else {
                    format!("{} -> {}", ga.name, gb.name)
                };
                entries.push(Entry {
                    name,
                    change: compare(&ga.atoms, &gb.atoms),
                    len_a: Some(ga.len()),
                    len_b: Some(gb.len()),
                });
            }
            (Some(ga), None) => entries.push(Entry {
                name: ga.name.clone(),
                change: Change::Removed,
                len_a: Some(ga.len()),
                len_b: None,
            }),
            (None, Some(gb)) => entries.push(Entry {
                name: gb.name.clone(),
                change: Change::Added,
                len_a: None,
                len_b: Some(gb.len()),
            }),
            (None, None) => unreachable!("i < max(len)"),
        }
    }
    Diff { entries }
}

fn diff_by_name(a: &IndexFile, b: &IndexFile) -> Diff {
    let mut entries = Vec::new();
    // Duplicate names are legal, so pair the k-th "SOL" in A with the k-th "SOL" in B.
    let mut used_b = vec![false; b.len()];

    for ga in &a.groups {
        let partner = b
            .groups
            .iter()
            .enumerate()
            .find(|(j, gb)| !used_b[*j] && gb.name == ga.name)
            .map(|(j, gb)| {
                used_b[j] = true;
                gb
            });
        match partner {
            Some(gb) => entries.push(Entry {
                name: ga.name.clone(),
                change: compare(&ga.atoms, &gb.atoms),
                len_a: Some(ga.len()),
                len_b: Some(gb.len()),
            }),
            None => entries.push(Entry {
                name: ga.name.clone(),
                change: Change::Removed,
                len_a: Some(ga.len()),
                len_b: None,
            }),
        }
    }

    for (j, gb) in b.groups.iter().enumerate() {
        if !used_b[j] {
            entries.push(Entry {
                name: gb.name.clone(),
                change: Change::Added,
                len_a: None,
                len_b: Some(gb.len()),
            });
        }
    }

    Diff { entries }
}

fn compare(a: &[u32], b: &[u32]) -> Change {
    if a == b {
        return Change::Same;
    }
    let (sa, sb) = (
        AtomSet::from_unsorted(a.to_vec()),
        AtomSet::from_unsorted(b.to_vec()),
    );
    if sa == sb {
        // Same atoms, different sequence. Worth flagging: index order matters to some gmx tools.
        return Change::Reordered;
    }
    Change::Changed {
        added: sb.difference(&sa),
        removed: sa.difference(&sb),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Group;

    fn f(groups: Vec<Group>) -> IndexFile {
        IndexFile { groups }
    }

    #[test]
    fn identical_files_are_empty() {
        let a = f(vec![Group::new("A", vec![1, 2])]);
        let d = diff(&a, &a, Match::Name);
        assert!(d.is_empty());
    }

    #[test]
    fn added_and_removed_groups() {
        let a = f(vec![Group::new("Gone", vec![1])]);
        let b = f(vec![Group::new("New", vec![2])]);
        let d = diff(&a, &b, Match::Name);
        assert_eq!(d.entries[0].change, Change::Removed);
        assert_eq!(d.entries[1].change, Change::Added);
        assert!(!d.is_empty());
    }

    #[test]
    fn changed_atoms_are_itemized() {
        let a = f(vec![Group::new("A", vec![1, 2, 3])]);
        let b = f(vec![Group::new("A", vec![2, 3, 4])]);
        let d = diff(&a, &b, Match::Name);
        match &d.entries[0].change {
            Change::Changed { added, removed } => {
                assert_eq!(added.as_slice(), [4]);
                assert_eq!(removed.as_slice(), [1]);
            }
            other => panic!("expected Changed, got {other:?}"),
        }
    }

    #[test]
    fn reordering_is_detected_not_hidden() {
        let a = f(vec![Group::new("A", vec![1, 2, 3])]);
        let b = f(vec![Group::new("A", vec![3, 2, 1])]);
        let d = diff(&a, &b, Match::Name);
        assert_eq!(d.entries[0].change, Change::Reordered);
        assert!(!d.is_empty());
    }

    #[test]
    fn duplicate_atoms_count_as_a_difference() {
        let a = f(vec![Group::new("A", vec![1, 1, 2])]);
        let b = f(vec![Group::new("A", vec![1, 2])]);
        let d = diff(&a, &b, Match::Name);
        assert_eq!(d.entries[0].change, Change::Reordered);
    }

    #[test]
    fn duplicate_names_pair_up_in_order() {
        let a = f(vec![
            Group::new("SOL", vec![1]),
            Group::new("SOL", vec![2]),
        ]);
        let b = f(vec![
            Group::new("SOL", vec![1]),
            Group::new("SOL", vec![9]),
        ]);
        let d = diff(&a, &b, Match::Name);
        assert_eq!(d.entries[0].change, Change::Same);
        assert!(matches!(d.entries[1].change, Change::Changed { .. }));
    }

    #[test]
    fn name_matching_ignores_position() {
        let a = f(vec![Group::new("A", vec![1]), Group::new("B", vec![2])]);
        let b = f(vec![Group::new("B", vec![2]), Group::new("A", vec![1])]);
        assert!(diff(&a, &b, Match::Name).is_empty());
    }

    #[test]
    fn index_matching_sees_the_swap() {
        let a = f(vec![Group::new("A", vec![1]), Group::new("B", vec![2])]);
        let b = f(vec![Group::new("B", vec![2]), Group::new("A", vec![1])]);
        let d = diff(&a, &b, Match::Index);
        assert!(!d.is_empty());
        assert_eq!(d.entries[0].name, "A -> B");
    }

    #[test]
    fn a_rename_shows_up_as_add_plus_remove_by_name() {
        let a = f(vec![Group::new("Old", vec![1])]);
        let b = f(vec![Group::new("New", vec![1])]);
        let d = diff(&a, &b, Match::Name);
        assert_eq!(d.entries.len(), 2);
    }
}
