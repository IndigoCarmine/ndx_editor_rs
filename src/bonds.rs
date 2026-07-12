//! The bond graph: what "adjacent" means.
//!
//! Two sources. A `.top` gives the *real* bonds, expanded across `[ molecules ]`. Without one, a
//! `.gro` alone can only estimate them from interatomic distances and covalent radii — good enough
//! to find the hydrogens hanging off a carbon, but a guess, and reported as such.

use std::collections::HashMap;

use crate::atomset::AtomSet;
use crate::model::AtomId;
use crate::structure::Structure;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondSource {
    /// `[ bonds ]` / `[ constraints ]` / `[ settles ]` out of a topology.
    Topology,
    /// Estimated from coordinates and covalent radii.
    Distance,
}

/// Adjacency in compressed-sparse-row form: one contiguous neighbour list, indexed by atom.
#[derive(Clone, Debug)]
pub struct BondGraph {
    offsets: Vec<u32>,
    neighbors: Vec<AtomId>,
    pub source: BondSource,
    pub nbonds: usize,
}

impl BondGraph {
    pub fn from_pairs(natoms: u32, pairs: &[(AtomId, AtomId)], source: BondSource) -> Self {
        let n = natoms as usize;
        let mut degree = vec![0u32; n + 1];
        let mut kept = 0usize;

        for &(a, b) in pairs {
            if a == 0 || b == 0 || a > natoms || b > natoms || a == b {
                continue;
            }
            degree[a as usize - 1] += 1;
            degree[b as usize - 1] += 1;
            kept += 1;
        }

        let mut offsets = vec![0u32; n + 1];
        let mut acc = 0u32;
        for i in 0..n {
            offsets[i] = acc;
            acc += degree[i];
        }
        offsets[n] = acc;

        let mut cursor = offsets.clone();
        let mut neighbors = vec![0 as AtomId; acc as usize];
        for &(a, b) in pairs {
            if a == 0 || b == 0 || a > natoms || b > natoms || a == b {
                continue;
            }
            neighbors[cursor[a as usize - 1] as usize] = b;
            cursor[a as usize - 1] += 1;
            neighbors[cursor[b as usize - 1] as usize] = a;
            cursor[b as usize - 1] += 1;
        }

        BondGraph {
            offsets,
            neighbors,
            source,
            nbonds: kept,
        }
    }

    pub fn natoms(&self) -> u32 {
        (self.offsets.len() - 1) as u32
    }

    pub fn neighbors(&self, a: AtomId) -> &[AtomId] {
        let Some(i) = (a as usize).checked_sub(1) else {
            return &[];
        };
        if i + 1 >= self.offsets.len() {
            return &[];
        }
        let (lo, hi) = (self.offsets[i] as usize, self.offsets[i + 1] as usize);
        &self.neighbors[lo..hi]
    }

    /// Every atom reachable from `set` in 1..=`depth` bonds.
    ///
    /// The result is the *neighbourhood*, so it contains an atom of `set` only when that atom is
    /// itself bonded to another atom of `set` — which is what makes `element H & bonded Protein`
    /// pick up the protein's own hydrogens. Write `bonded X & !X` for the strictly-outside version.
    pub fn expand(&self, set: &AtomSet, depth: u8) -> AtomSet {
        let mut frontier: Vec<AtomId> = set.as_slice().to_vec();
        let mut reached: Vec<AtomId> = Vec::new();
        let mut seen: Vec<bool> = vec![false; self.natoms() as usize + 1];

        for _ in 0..depth.max(1) {
            let mut next: Vec<AtomId> = Vec::new();
            for &a in &frontier {
                for &b in self.neighbors(a) {
                    if !seen[b as usize] {
                        seen[b as usize] = true;
                        reached.push(b);
                        next.push(b);
                    }
                }
            }
            if next.is_empty() {
                break;
            }
            frontier = next;
        }

        AtomSet::from_unsorted(reached)
    }
}

/// Covalent radii in nm. Anything unlisted gets carbon's, which errs towards finding a bond
/// rather than missing one.
fn covalent_radius_nm(element: &str) -> f32 {
    match element {
        "H" => 0.031,
        "C" => 0.076,
        "N" => 0.071,
        "O" => 0.066,
        "F" => 0.057,
        "P" => 0.107,
        "S" => 0.105,
        "CL" => 0.102,
        "BR" => 0.120,
        "I" => 0.139,
        "NA" => 0.166,
        "MG" => 0.141,
        "K" => 0.203,
        "CA" => 0.176,
        "ZN" => 0.122,
        "FE" => 0.132,
        "CU" => 0.132,
        _ => 0.077,
    }
}

const EXTRA_TOLERANCE: f32 = 0.045;
const MIN_DISTANCE: f32 = 0.04;
const MAX_RADIUS: f32 = 0.203;

/// Estimate bonds from coordinates, using a uniform grid so it stays linear in the atom count.
///
/// Two atoms are bonded when they are closer than `r_i + r_j + 0.045 nm` (and not on top of each
/// other, which would be a duplicate atom rather than a bond).
pub fn infer_from_distance(s: &Structure) -> BondGraph {
    let cell_size = MAX_RADIUS * 2.0 + EXTRA_TOLERANCE;
    let inv = 1.0 / cell_size;
    let min_d2 = MIN_DISTANCE * MIN_DISTANCE;

    let mut grid: HashMap<(i32, i32, i32), Vec<u32>> = HashMap::new();
    let mut pairs: Vec<(AtomId, AtomId)> = Vec::new();

    for i in 0..s.positions.len() {
        let p = s.positions[i];
        let cell = (
            (p[0] * inv).floor() as i32,
            (p[1] * inv).floor() as i32,
            (p[2] * inv).floor() as i32,
        );
        let ri = covalent_radius_nm(&s.elements[i]);

        // Only the 27 cells around this atom can hold a partner, since a bond is shorter than one
        // cell by construction.
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(bucket) = grid.get(&(cell.0 + dx, cell.1 + dy, cell.2 + dz)) else {
                        continue;
                    };
                    for &j in bucket {
                        let q = s.positions[j as usize];
                        let d2 = (q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2) + (q[2] - p[2]).powi(2);
                        if d2 < min_d2 {
                            continue;
                        }
                        let rj = covalent_radius_nm(&s.elements[j as usize]);
                        let max = ri + rj + EXTRA_TOLERANCE;
                        if d2 <= max * max {
                            // 1-based, and the bucket only ever holds earlier atoms, so each pair
                            // is produced exactly once.
                            pairs.push((j + 1, i as u32 + 1));
                        }
                    }
                }
            }
        }

        grid.entry(cell).or_default().push(i as u32);
    }

    BondGraph::from_pairs(s.natoms(), &pairs, BondSource::Distance)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Ethane-ish: C1-C2, each carbon carrying three hydrogens.
    fn ethane_pairs() -> Vec<(AtomId, AtomId)> {
        vec![
            (1, 2),
            (1, 3),
            (1, 4),
            (1, 5), // C1 - H,H,H and C2
            (5, 6),
            (5, 7),
            (5, 8), // C2 - H,H,H
        ]
    }

    fn graph() -> BondGraph {
        BondGraph::from_pairs(8, &ethane_pairs(), BondSource::Topology)
    }

    #[test]
    fn neighbours_are_symmetric() {
        let g = graph();
        assert!(g.neighbors(1).contains(&5));
        assert!(g.neighbors(5).contains(&1));
        assert_eq!(g.neighbors(2), [1]);
        assert_eq!(g.nbonds, 7);
    }

    #[test]
    fn out_of_range_atoms_have_no_neighbours() {
        let g = graph();
        assert!(g.neighbors(0).is_empty());
        assert!(g.neighbors(99).is_empty());
    }

    #[test]
    fn bad_pairs_are_dropped_rather_than_panicking() {
        let g = BondGraph::from_pairs(3, &[(1, 2), (0, 1), (1, 99), (2, 2)], BondSource::Topology);
        assert_eq!(g.nbonds, 1);
    }

    #[test]
    fn expand_finds_the_neighbourhood() {
        let g = graph();
        // The atoms bonded to C1 = {2,3,4,5}: its three hydrogens plus the other carbon.
        let set = AtomSet::from_unsorted(vec![1]);
        assert_eq!(g.expand(&set, 1).as_slice(), [2, 3, 4, 5]);
    }

    /// This is the shape of the target query: the H atoms attached to a group of carbons.
    #[test]
    fn hydrogens_attached_to_a_carbon() {
        let g = graph();
        let carbons = AtomSet::from_unsorted(vec![1, 5]);
        let neighbourhood = g.expand(&carbons, 1);
        // Every hydrogen, plus each carbon (they are bonded to each other).
        assert_eq!(neighbourhood.as_slice(), [1, 2, 3, 4, 5, 6, 7, 8]);
        // `element H & bonded C` then narrows it; `bonded X & !X` is the strictly-outside form.
        let outside = neighbourhood.difference(&carbons);
        assert_eq!(outside.as_slice(), [2, 3, 4, 6, 7, 8]);
    }

    #[test]
    fn depth_walks_further() {
        let g = graph();
        let h = AtomSet::from_unsorted(vec![2]); // one hydrogen on C1
        assert_eq!(g.expand(&h, 1).as_slice(), [1]);
        // Two bonds out: C1, then C1's other neighbours.
        assert_eq!(g.expand(&h, 2).as_slice(), [1, 2, 3, 4, 5]);
        assert_eq!(g.expand(&h, 3).as_slice(), [1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn expanding_an_empty_set_yields_nothing() {
        assert!(graph().expand(&AtomSet::empty(), 2).is_empty());
    }

    #[test]
    fn distance_inference_finds_a_water() {
        // O at the origin, two H at 0.1 nm — a real O-H bond is 0.096 nm.
        let s = Structure {
            names: vec!["OW".into(), "HW1".into(), "HW2".into()],
            resnames: vec!["SOL".into(); 3],
            resids: vec![1; 3],
            elements: vec!["O".into(), "H".into(), "H".into()],
            positions: vec![[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.1, 0.0]],
            ..Default::default()
        };
        let g = infer_from_distance(&s);
        assert_eq!(g.source, BondSource::Distance);
        assert_eq!(g.nbonds, 2);
        assert_eq!(g.neighbors(1), [2, 3]);
        // The two hydrogens are 0.141 nm apart — far past H+H+tol = 0.107 — so no H-H bond.
        assert_eq!(g.neighbors(2), [1]);
    }

    #[test]
    fn distance_inference_ignores_distant_atoms() {
        let s = Structure {
            names: vec!["C".into(), "C".into()],
            resnames: vec!["X".into(); 2],
            resids: vec![1; 2],
            elements: vec!["C".into(), "C".into()],
            positions: vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]],
            ..Default::default()
        };
        assert_eq!(infer_from_distance(&s).nbonds, 0);
    }

    /// Atoms sitting on top of each other are a duplicated coordinate, not a bond.
    #[test]
    fn distance_inference_skips_coincident_atoms() {
        let s = Structure {
            names: vec!["C".into(), "C".into()],
            resnames: vec!["X".into(); 2],
            resids: vec![1; 2],
            elements: vec!["C".into(), "C".into()],
            positions: vec![[0.0; 3], [0.0; 3]],
            ..Default::default()
        };
        assert_eq!(infer_from_distance(&s).nbonds, 0);
    }
}
