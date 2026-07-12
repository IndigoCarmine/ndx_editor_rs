use std::fmt;
use std::str::FromStr;

use crate::atomset::AtomSet;
use crate::error::{NdxError, Result};

/// Atom index as it appears in the file: 1-based.
pub type AtomId = u32;

/// One `[ name ]` block.
///
/// `atoms` is kept in **file order, with duplicates**, so that reading and writing a file
/// round-trips byte-for-byte. Set operations go through [`AtomSet`], which normalizes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Group {
    pub name: String,
    pub atoms: Vec<AtomId>,
}

impl Group {
    pub fn new(name: impl Into<String>, atoms: Vec<AtomId>) -> Self {
        Group {
            name: name.into(),
            atoms,
        }
    }

    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    pub fn is_sorted(&self) -> bool {
        self.atoms.windows(2).all(|w| w[0] <= w[1])
    }

    pub fn has_duplicates(&self) -> bool {
        let mut sorted = self.atoms.clone();
        sorted.sort_unstable();
        let n = sorted.len();
        sorted.dedup();
        sorted.len() != n
    }

    /// Sorted, deduplicated view for set arithmetic.
    pub fn to_set(&self) -> AtomSet {
        AtomSet::from_unsorted(self.atoms.clone())
    }

    pub fn max_atom(&self) -> Option<AtomId> {
        self.atoms.iter().copied().max()
    }
}

/// A whole `.ndx` file. A group's position in `groups` is the 0-based id used in expressions.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct IndexFile {
    pub groups: Vec<Group>,
}

impl IndexFile {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.groups.len()
    }

    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }

    pub fn get(&self, id: usize) -> Result<&Group> {
        self.groups.get(id).ok_or(NdxError::GroupIdOutOfRange {
            id,
            count: self.groups.len(),
        })
    }

    pub fn get_mut(&mut self, id: usize) -> Result<&mut Group> {
        let count = self.groups.len();
        self.groups
            .get_mut(id)
            .ok_or(NdxError::GroupIdOutOfRange { id, count })
    }

    /// Duplicate group names are legal in `.ndx`, so this returns every match.
    pub fn find_by_name(&self, name: &str) -> Vec<usize> {
        self.groups
            .iter()
            .enumerate()
            .filter(|(_, g)| g.name == name)
            .map(|(i, _)| i)
            .collect()
    }

    /// Every group a reference names. A glob may name several; anything else names exactly one.
    ///
    /// Results are in file order, and never empty — a pattern that matches nothing is an error,
    /// not a silent no-op.
    pub fn resolve_all(&self, r: &GroupRef) -> Result<Vec<usize>> {
        match r {
            GroupRef::Id(id) => {
                if *id < self.groups.len() {
                    Ok(vec![*id])
                } else {
                    Err(NdxError::GroupIdOutOfRange {
                        id: *id,
                        count: self.groups.len(),
                    })
                }
            }

            GroupRef::Glob(pattern) => {
                let hits: Vec<usize> = self
                    .groups
                    .iter()
                    .enumerate()
                    .filter(|(_, g)| crate::glob::matches(pattern, &g.name))
                    .map(|(i, _)| i)
                    .collect();
                if hits.is_empty() {
                    return Err(NdxError::NoGroupMatches {
                        pattern: pattern.clone(),
                    });
                }
                Ok(hits)
            }

            GroupRef::Name { name, occurrence } => {
                let hits = self.find_by_name(name);
                match (hits.len(), occurrence) {
                    (0, _) => Err(NdxError::UnknownGroup {
                        name: name.clone(),
                        suggestion: self.suggest_name(name).map(str::to_owned),
                    }),
                    (_, Some(k)) => hits
                        .get(k.saturating_sub(1))
                        .copied()
                        .map(|i| vec![i])
                        .ok_or_else(|| NdxError::NoSuchOccurrence {
                            name: name.clone(),
                            have: hits.len(),
                            want: *k,
                        }),
                    (1, None) => Ok(vec![hits[0]]),
                    (_, None) => Err(NdxError::AmbiguousName {
                        name: name.clone(),
                        ids: hits,
                    }),
                }
            }
        }
    }

    /// The single group a reference names. Errors when a glob matched more than one, because the
    /// caller (rename, split, --universe) can only act on one.
    pub fn resolve(&self, r: &GroupRef) -> Result<usize> {
        let hits = self.resolve_all(r)?;
        match hits.len() {
            1 => Ok(hits[0]),
            _ => Err(NdxError::PatternNotUnique {
                pattern: r.to_string(),
                names: hits.iter().map(|&i| self.groups[i].name.clone()).collect(),
            }),
        }
    }

    /// Resolve several references, flattening globs. Duplicates are removed, file order is kept.
    pub fn resolve_many(&self, refs: &[GroupRef]) -> Result<Vec<usize>> {
        let mut out: Vec<usize> = Vec::new();
        for r in refs {
            for id in self.resolve_all(r)? {
                if !out.contains(&id) {
                    out.push(id);
                }
            }
        }
        Ok(out)
    }

    pub fn push(&mut self, g: Group) -> usize {
        self.groups.push(g);
        self.groups.len() - 1
    }

    /// Remove the given ids. Remaining groups keep their relative order and are renumbered.
    pub fn remove_many(&mut self, ids: &[usize]) {
        let mut doomed: Vec<usize> = ids.to_vec();
        doomed.sort_unstable();
        doomed.dedup();
        for id in doomed.into_iter().rev() {
            if id < self.groups.len() {
                self.groups.remove(id);
            }
        }
    }

    /// Keep only the given ids, **in the order given** (so `keep` doubles as a reorder).
    pub fn retain_ids(&mut self, ids: &[usize]) {
        let kept: Vec<Group> = ids
            .iter()
            .filter_map(|&i| self.groups.get(i).cloned())
            .collect();
        self.groups = kept;
    }

    pub fn union_all(&self) -> AtomSet {
        let mut all: Vec<AtomId> = Vec::new();
        for g in &self.groups {
            all.extend_from_slice(&g.atoms);
        }
        AtomSet::from_unsorted(all)
    }

    pub fn max_atom(&self) -> Option<AtomId> {
        self.groups.iter().filter_map(Group::max_atom).max()
    }

    /// Closest group name by edit distance, for "did you mean" hints.
    pub fn suggest_name(&self, miss: &str) -> Option<&str> {
        // Tolerate more typos in longer names, but never suggest something wildly different.
        let budget = (miss.chars().count() / 3).clamp(1, 3);
        self.groups
            .iter()
            .map(|g| (edit_distance(&g.name.to_lowercase(), &miss.to_lowercase()), g))
            .filter(|(d, _)| *d <= budget)
            .min_by_key(|(d, _)| *d)
            .map(|(_, g)| g.name.as_str())
    }

    /// Make `desired` unique among existing group names by appending `_2`, `_3`, ...
    pub fn unique_name(&self, desired: &str) -> String {
        if self.find_by_name(desired).is_empty() {
            return desired.to_string();
        }
        (2..)
            .map(|n| format!("{desired}_{n}"))
            .find(|cand| self.find_by_name(cand).is_empty())
            .expect("infinite range yields a free name")
    }
}

/// Levenshtein distance, two-row variant.
fn edit_distance(a: &str, b: &str) -> usize {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut cur = vec![0usize; b.len() + 1];
    for (i, ca) in a.iter().enumerate() {
        cur[0] = i + 1;
        for (j, cb) in b.iter().enumerate() {
            let cost = usize::from(ca != cb);
            cur[j + 1] = (prev[j + 1] + 1).min(cur[j] + 1).min(prev[j] + cost);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[b.len()]
}

/// How the user names a group: by 0-based id, by name (with `#K` to pick among duplicates), or by
/// a glob pattern that may match several.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GroupRef {
    Id(usize),
    Name {
        name: String,
        occurrence: Option<usize>,
    },
    /// `Fiber*`, `O?`. Matches every group whose name fits, in file order.
    Glob(String),
}

impl GroupRef {
    pub fn name(name: impl Into<String>) -> Self {
        GroupRef::Name {
            name: name.into(),
            occurrence: None,
        }
    }

    /// A bare word from the user: a glob if it has `*`/`?` in it, an exact name otherwise.
    pub fn from_bare(s: impl Into<String>) -> Self {
        let s = s.into();
        if crate::glob::is_pattern(&s) {
            GroupRef::Glob(s)
        } else {
            GroupRef::name(s)
        }
    }
}

impl FromStr for GroupRef {
    type Err = NdxError;

    /// `3` | `SOL` | `SOL#2` | `Fiber*` | `"C alpha"` | `"literal*star"`
    ///
    /// Quoting turns off both the id and the glob reading, so a group genuinely called `*` or `13`
    /// is still reachable.
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        let s = s.trim();
        if s.is_empty() {
            return Err(NdxError::Other("empty group reference".into()));
        }
        if let Some(inner) = unquote(s) {
            return Ok(GroupRef::name(inner));
        }
        if let Ok(id) = s.parse::<usize>() {
            return Ok(GroupRef::Id(id));
        }
        // Split a trailing `#K`, but only when K is a number — group names may contain '#'.
        if let Some((head, tail)) = s.rsplit_once('#')
            && let Ok(k) = tail.parse::<usize>()
            && k > 0
        {
            return Ok(GroupRef::Name {
                name: unquote(head).unwrap_or(head).to_string(),
                occurrence: Some(k),
            });
        }
        Ok(GroupRef::from_bare(s))
    }
}

/// `Some(inner)` when `s` was wrapped in double quotes.
fn unquote(s: &str) -> Option<&str> {
    s.strip_prefix('"').and_then(|s| s.strip_suffix('"'))
}

impl fmt::Display for GroupRef {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GroupRef::Id(id) => write!(f, "{id}"),
            GroupRef::Glob(p) => write!(f, "{p}"),
            GroupRef::Name {
                name,
                occurrence: None,
            } => write!(f, "{name}"),
            GroupRef::Name {
                name,
                occurrence: Some(k),
            } => write!(f, "{name}#{k}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> IndexFile {
        IndexFile {
            groups: vec![
                Group::new("System", vec![1, 2, 3, 4, 5]),
                Group::new("Protein", vec![1, 2]),
                Group::new("SOL", vec![3, 4, 5]),
                Group::new("SOL", vec![3, 4]),
            ],
        }
    }

    #[test]
    fn resolve_by_id() {
        assert_eq!(sample().resolve(&GroupRef::Id(1)).unwrap(), 1);
    }

    #[test]
    fn resolve_out_of_range() {
        assert!(matches!(
            sample().resolve(&GroupRef::Id(9)),
            Err(NdxError::GroupIdOutOfRange { id: 9, count: 4 })
        ));
    }

    #[test]
    fn resolve_unique_name() {
        assert_eq!(sample().resolve(&GroupRef::name("Protein")).unwrap(), 1);
    }

    #[test]
    fn duplicate_name_is_ambiguous() {
        assert!(matches!(
            sample().resolve(&GroupRef::name("SOL")),
            Err(NdxError::AmbiguousName { ref ids, .. }) if ids == &[2, 3]
        ));
    }

    #[test]
    fn occurrence_disambiguates() {
        let r = GroupRef::Name {
            name: "SOL".into(),
            occurrence: Some(2),
        };
        assert_eq!(sample().resolve(&r).unwrap(), 3);
    }

    #[test]
    fn occurrence_out_of_range() {
        let r = GroupRef::Name {
            name: "SOL".into(),
            occurrence: Some(3),
        };
        assert!(matches!(
            sample().resolve(&r),
            Err(NdxError::NoSuchOccurrence { have: 2, want: 3, .. })
        ));
    }

    #[test]
    fn unknown_name_suggests() {
        let err = sample().resolve(&GroupRef::name("Protien")).unwrap_err();
        match err {
            NdxError::UnknownGroup { suggestion, .. } => {
                assert_eq!(suggestion.as_deref(), Some("Protein"));
            }
            other => panic!("expected UnknownGroup, got {other:?}"),
        }
    }

    #[test]
    fn wildly_different_name_suggests_nothing() {
        let err = sample().resolve(&GroupRef::name("zzzzzzzz")).unwrap_err();
        assert!(matches!(
            err,
            NdxError::UnknownGroup {
                suggestion: None,
                ..
            }
        ));
    }

    fn globby() -> IndexFile {
        IndexFile {
            groups: vec![
                Group::new("Fiber1", vec![1, 2]),
                Group::new("Fiber2", vec![3]),
                Group::new("fiberA", vec![4]),
                Group::new("AlkylA", vec![5]),
                Group::new("OA", vec![6]),
                Group::new("OB", vec![7]),
            ],
        }
    }

    #[test]
    fn a_glob_resolves_to_every_match_in_file_order() {
        let f = globby();
        let r: GroupRef = "Fiber*".parse().unwrap();
        assert_eq!(f.resolve_all(&r).unwrap(), [0, 1]);

        let r: GroupRef = "O*".parse().unwrap();
        assert_eq!(f.resolve_all(&r).unwrap(), [4, 5]);

        // Case-sensitive: `Fiber*` must not pick up `fiberA`.
        let r: GroupRef = "*A".parse().unwrap();
        assert_eq!(f.resolve_all(&r).unwrap(), [2, 3, 4]);
    }

    #[test]
    fn a_glob_matching_nothing_is_an_error_not_a_silent_no_op() {
        assert!(matches!(
            globby().resolve_all(&"Nope*".parse().unwrap()),
            Err(NdxError::NoGroupMatches { .. })
        ));
    }

    #[test]
    fn resolve_needs_a_glob_to_be_unique() {
        let f = globby();
        // rename/split/--universe can only act on one group.
        assert!(matches!(
            f.resolve(&"Fiber*".parse().unwrap()),
            Err(NdxError::PatternNotUnique { ref names, .. }) if names == &["Fiber1", "Fiber2"]
        ));
        // A glob with exactly one match is fine.
        assert_eq!(f.resolve(&"Alkyl*".parse().unwrap()).unwrap(), 3);
    }

    #[test]
    fn resolve_many_flattens_globs_and_dedups() {
        let f = globby();
        let refs: Vec<GroupRef> = ["O*", "Fiber1", "OA"]
            .iter()
            .map(|s| s.parse().unwrap())
            .collect();
        // OA appears twice across the refs but must be listed once, and order-of-mention wins.
        assert_eq!(f.resolve_many(&refs).unwrap(), [4, 5, 0]);
    }

    #[test]
    fn quoting_makes_a_pattern_literal() {
        let f = IndexFile {
            groups: vec![Group::new("O*", vec![1]), Group::new("OA", vec![2])],
        };
        // Bare: a pattern, matching both.
        assert_eq!(f.resolve_all(&"O*".parse().unwrap()).unwrap(), [0, 1]);
        // Quoted: the group actually named "O*".
        assert_eq!(f.resolve_all(&"\"O*\"".parse().unwrap()).unwrap(), [0]);
    }

    #[test]
    fn group_ref_parsing() {
        assert_eq!("3".parse::<GroupRef>().unwrap(), GroupRef::Id(3));
        assert_eq!(
            "SOL".parse::<GroupRef>().unwrap(),
            GroupRef::name("SOL")
        );
        assert_eq!(
            "SOL#2".parse::<GroupRef>().unwrap(),
            GroupRef::Name {
                name: "SOL".into(),
                occurrence: Some(2)
            }
        );
        assert_eq!(
            "\"C alpha\"".parse::<GroupRef>().unwrap(),
            GroupRef::name("C alpha")
        );
        // A '#' that isn't followed by a number is part of the name.
        assert_eq!("a#b".parse::<GroupRef>().unwrap(), GroupRef::name("a#b"));
    }

    #[test]
    fn remove_many_renumbers() {
        let mut f = sample();
        f.remove_many(&[1, 3]);
        let names: Vec<&str> = f.groups.iter().map(|g| g.name.as_str()).collect();
        assert_eq!(names, ["System", "SOL"]);
    }

    #[test]
    fn retain_ids_reorders() {
        let mut f = sample();
        f.retain_ids(&[2, 0]);
        let names: Vec<&str> = f.groups.iter().map(|g| g.name.as_str()).collect();
        assert_eq!(names, ["SOL", "System"]);
    }

    #[test]
    fn unique_name_appends_suffix() {
        let f = sample();
        assert_eq!(f.unique_name("Protein"), "Protein_2");
        assert_eq!(f.unique_name("Fresh"), "Fresh");
    }

    #[test]
    fn flags_detect_unsorted_and_dups() {
        let g = Group::new("x", vec![3, 1, 1]);
        assert!(!g.is_sorted());
        assert!(g.has_duplicates());
        assert_eq!(g.to_set().as_slice(), [1, 3]);
    }
}
