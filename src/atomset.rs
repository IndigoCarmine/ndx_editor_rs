use crate::error::{NdxError, Result};
use crate::model::AtomId;

/// A set of atoms. Invariant: strictly ascending (sorted, no duplicates).
///
/// This is the *computation* type. Groups keep their file order in a plain `Vec`; only set
/// arithmetic goes through here, and it always yields a normalized result — same as make_ndx.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AtomSet(Vec<AtomId>);

impl AtomSet {
    pub fn empty() -> Self {
        AtomSet(Vec::new())
    }

    pub fn from_unsorted(mut v: Vec<AtomId>) -> Self {
        v.sort_unstable();
        v.dedup();
        AtomSet(v)
    }

    pub fn from_sorted_unique(v: Vec<AtomId>) -> Self {
        debug_assert!(
            v.windows(2).all(|w| w[0] < w[1]),
            "AtomSet::from_sorted_unique got a non-ascending slice"
        );
        AtomSet(v)
    }

    /// `lo..=hi`, empty when `lo > hi`.
    pub fn range_inclusive(lo: AtomId, hi: AtomId) -> Self {
        AtomSet(if lo > hi { Vec::new() } else { (lo..=hi).collect() })
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn contains(&self, a: AtomId) -> bool {
        self.0.binary_search(&a).is_ok()
    }

    pub fn as_slice(&self) -> &[AtomId] {
        &self.0
    }

    pub fn into_vec(self) -> Vec<AtomId> {
        self.0
    }

    pub fn max(&self) -> Option<AtomId> {
        self.0.last().copied()
    }

    pub fn union(&self, other: &Self) -> Self {
        let (a, b) = (&self.0, &other.0);
        let mut out = Vec::with_capacity(a.len() + b.len());
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            match a[i].cmp(&b[j]) {
                std::cmp::Ordering::Less => {
                    out.push(a[i]);
                    i += 1;
                }
                std::cmp::Ordering::Greater => {
                    out.push(b[j]);
                    j += 1;
                }
                std::cmp::Ordering::Equal => {
                    out.push(a[i]);
                    i += 1;
                    j += 1;
                }
            }
        }
        out.extend_from_slice(&a[i..]);
        out.extend_from_slice(&b[j..]);
        AtomSet(out)
    }

    pub fn intersection(&self, other: &Self) -> Self {
        let (a, b) = (&self.0, &other.0);
        let mut out = Vec::with_capacity(a.len().min(b.len()));
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            match a[i].cmp(&b[j]) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    out.push(a[i]);
                    i += 1;
                    j += 1;
                }
            }
        }
        AtomSet(out)
    }

    pub fn difference(&self, other: &Self) -> Self {
        let (a, b) = (&self.0, &other.0);
        let mut out = Vec::with_capacity(a.len());
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            match a[i].cmp(&b[j]) {
                std::cmp::Ordering::Less => {
                    out.push(a[i]);
                    i += 1;
                }
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    i += 1;
                    j += 1;
                }
            }
        }
        out.extend_from_slice(&a[i..]);
        AtomSet(out)
    }

    pub fn complement(&self, universe: &Self) -> Self {
        universe.difference(self)
    }

    pub fn is_subset_of(&self, other: &Self) -> bool {
        self.difference(other).is_empty()
    }
}

impl FromIterator<AtomId> for AtomSet {
    fn from_iter<I: IntoIterator<Item = AtomId>>(iter: I) -> Self {
        AtomSet::from_unsorted(iter.into_iter().collect())
    }
}

/// Compress a set into inclusive ranges: `[1,2,3,7,9,10]` -> `[(1,3), (7,7), (9,10)]`.
pub fn to_ranges(s: &AtomSet) -> Vec<(AtomId, AtomId)> {
    let mut out: Vec<(AtomId, AtomId)> = Vec::new();
    for &a in s.as_slice() {
        match out.last_mut() {
            Some(last) if last.1 + 1 == a => last.1 = a,
            _ => out.push((a, a)),
        }
    }
    out
}

/// `"1-10,15,20-30"`
pub fn format_ranges(s: &AtomSet) -> String {
    to_ranges(s)
        .into_iter()
        .map(|(lo, hi)| {
            if lo == hi {
                lo.to_string()
            } else {
                format!("{lo}-{hi}")
            }
        })
        .collect::<Vec<_>>()
        .join(",")
}

/// Parse `"1-10,15,20-30"` (also accepts whitespace as a separator).
pub fn parse_ranges(s: &str) -> Result<AtomSet> {
    let mut atoms: Vec<AtomId> = Vec::new();
    for tok in s.split([',', ' ', '\t']).filter(|t| !t.is_empty()) {
        let (lo, hi) = match tok.split_once('-') {
            Some((lo, hi)) => (parse_atom(lo)?, parse_atom(hi)?),
            None => {
                let a = parse_atom(tok)?;
                (a, a)
            }
        };
        if lo > hi {
            return Err(NdxError::Other(format!(
                "range {tok:?} runs backwards ({lo} > {hi})"
            )));
        }
        atoms.extend(lo..=hi);
    }
    Ok(AtomSet::from_unsorted(atoms))
}

fn parse_atom(t: &str) -> Result<AtomId> {
    let t = t.trim();
    let a: AtomId = t
        .parse()
        .map_err(|_| NdxError::Other(format!("{t:?} is not an atom index")))?;
    if a == 0 {
        return Err(NdxError::Other(
            "atom index 0 is invalid (.ndx indices are 1-based)".into(),
        ));
    }
    Ok(a)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set(v: &[AtomId]) -> AtomSet {
        AtomSet::from_unsorted(v.to_vec())
    }

    #[test]
    fn from_unsorted_normalizes() {
        assert_eq!(set(&[3, 1, 2, 1]).as_slice(), [1, 2, 3]);
    }

    #[test]
    fn union_merges() {
        assert_eq!(set(&[1, 3, 5]).union(&set(&[2, 3, 6])).as_slice(), [1, 2, 3, 5, 6]);
    }

    #[test]
    fn union_with_empty() {
        assert_eq!(set(&[1, 2]).union(&AtomSet::empty()).as_slice(), [1, 2]);
        assert_eq!(AtomSet::empty().union(&set(&[1, 2])).as_slice(), [1, 2]);
    }

    #[test]
    fn intersection_of_disjoint_is_empty() {
        assert!(set(&[1, 3]).intersection(&set(&[2, 4])).is_empty());
    }

    #[test]
    fn intersection_of_overlap() {
        assert_eq!(set(&[1, 2, 3, 4]).intersection(&set(&[3, 4, 5])).as_slice(), [3, 4]);
    }

    #[test]
    fn difference_removes() {
        assert_eq!(set(&[1, 2, 3, 4]).difference(&set(&[2, 4])).as_slice(), [1, 3]);
    }

    #[test]
    fn difference_of_identical_is_empty() {
        assert!(set(&[1, 2, 3]).difference(&set(&[1, 2, 3])).is_empty());
    }

    #[test]
    fn complement_within_universe() {
        let u = AtomSet::range_inclusive(1, 5);
        assert_eq!(set(&[2, 4]).complement(&u).as_slice(), [1, 3, 5]);
    }

    #[test]
    fn range_inclusive_backwards_is_empty() {
        assert!(AtomSet::range_inclusive(5, 1).is_empty());
    }

    #[test]
    fn subset() {
        assert!(set(&[2, 3]).is_subset_of(&set(&[1, 2, 3, 4])));
        assert!(!set(&[2, 9]).is_subset_of(&set(&[1, 2, 3, 4])));
    }

    #[test]
    fn ranges_roundtrip() {
        let s = set(&[1, 2, 3, 7, 9, 10]);
        assert_eq!(format_ranges(&s), "1-3,7,9-10");
        assert_eq!(parse_ranges("1-3,7,9-10").unwrap(), s);
    }

    #[test]
    fn parse_ranges_accepts_whitespace_and_dedups() {
        assert_eq!(parse_ranges("1 2 2 3").unwrap().as_slice(), [1, 2, 3]);
    }

    #[test]
    fn parse_ranges_rejects_zero_and_backwards() {
        assert!(parse_ranges("0").is_err());
        assert!(parse_ranges("5-1").is_err());
        assert!(parse_ranges("a").is_err());
    }

    #[test]
    fn empty_set_formats_empty() {
        assert_eq!(format_ranges(&AtomSet::empty()), "");
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    fn any_set() -> impl Strategy<Value = AtomSet> {
        prop::collection::vec(1u32..60, 0..40).prop_map(AtomSet::from_unsorted)
    }

    fn ascending(s: &AtomSet) -> bool {
        s.as_slice().windows(2).all(|w| w[0] < w[1])
    }

    proptest! {
        #[test]
        fn results_are_always_strictly_ascending(a in any_set(), b in any_set()) {
            prop_assert!(ascending(&a.union(&b)));
            prop_assert!(ascending(&a.intersection(&b)));
            prop_assert!(ascending(&a.difference(&b)));
        }

        #[test]
        fn union_is_commutative(a in any_set(), b in any_set()) {
            prop_assert_eq!(a.union(&b), b.union(&a));
        }

        #[test]
        fn intersection_is_a_subset(a in any_set(), b in any_set()) {
            prop_assert!(a.intersection(&b).is_subset_of(&a));
        }

        #[test]
        fn inclusion_exclusion(a in any_set(), b in any_set()) {
            prop_assert_eq!(
                a.union(&b).len() + a.intersection(&b).len(),
                a.len() + b.len()
            );
        }

        /// The rewrite that keeps the universe out of the common path:
        /// for any universe U containing A and B, `A & !B == A \ B`.
        #[test]
        fn difference_equals_and_not(a in any_set(), b in any_set()) {
            let universe = AtomSet::range_inclusive(1, 60);
            prop_assert_eq!(a.intersection(&b.complement(&universe)), a.difference(&b));
        }

        #[test]
        fn double_complement_is_identity(a in any_set()) {
            let universe = AtomSet::range_inclusive(1, 60);
            prop_assert_eq!(a.complement(&universe).complement(&universe), a);
        }

        #[test]
        fn de_morgan(a in any_set(), b in any_set()) {
            let u = AtomSet::range_inclusive(1, 60);
            prop_assert_eq!(
                a.union(&b).complement(&u),
                a.complement(&u).intersection(&b.complement(&u))
            );
        }

        #[test]
        fn range_format_roundtrips(a in any_set()) {
            prop_assert_eq!(parse_ranges(&format_ranges(&a)).unwrap(), a);
        }
    }
}
