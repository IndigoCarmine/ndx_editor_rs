//! `within R of X`: every atom whose distance to the selection is at most R.

use std::collections::HashMap;

use crate::atomset::AtomSet;
use crate::model::AtomId;
use crate::structure::Structure;

/// Atoms within `radius` nm of any atom in `of`.
///
/// The selection is included in its own neighbourhood, since its atoms are at distance zero —
/// `within R of X & !X` is the exclusive version. With `pbc`, distances use the minimum-image
/// convention against the box from the `.gro`.
pub fn within(s: &Structure, of: &AtomSet, radius: f32, pbc: bool) -> AtomSet {
    if of.is_empty() || radius < 0.0 {
        return AtomSet::empty();
    }

    let cell = radius.max(1e-3);
    let inv = 1.0 / cell;
    let r2 = radius * radius;
    let box_diag = s.box_diag;
    let use_pbc = pbc && box_diag.iter().all(|v| *v > 1e-6);

    // Under PBC the grid has to wrap too, or a pair straddling the boundary lands in cells that
    // are numerically far apart and never gets compared — the distance check alone is not enough.
    // `ncells` is the number of cells per axis; the neighbour offsets are taken modulo it.
    let ncells: [i32; 3] = if use_pbc {
        [0, 1, 2].map(|k| ((box_diag[k] * inv).floor() as i32).max(1))
    } else {
        [0; 3]
    };
    let key = |p: [f32; 3]| cell_of(p, inv, use_pbc, box_diag, ncells);

    // Bucket the *selection*, then sweep every atom against it: the selection is usually the
    // smaller side, and this keeps one pass over the coordinates.
    let mut grid: HashMap<(i32, i32, i32), Vec<AtomId>> = HashMap::new();
    for &a in of.as_slice() {
        let Some(p) = s.pos(a) else { continue };
        grid.entry(key(p)).or_default().push(a);
    }

    let mut hits: Vec<AtomId> = Vec::new();
    for i in 0..s.positions.len() {
        let p = s.positions[i];
        let c = key(p);

        'atom: for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let nb = if use_pbc {
                        (
                            (c.0 + dx).rem_euclid(ncells[0]),
                            (c.1 + dy).rem_euclid(ncells[1]),
                            (c.2 + dz).rem_euclid(ncells[2]),
                        )
                    } else {
                        (c.0 + dx, c.1 + dy, c.2 + dz)
                    };
                    let Some(bucket) = grid.get(&nb) else { continue };
                    for &b in bucket {
                        let q = s.positions[b as usize - 1];
                        if dist2(p, q, box_diag, use_pbc) <= r2 {
                            hits.push(i as AtomId + 1);
                            break 'atom;
                        }
                    }
                }
            }
        }
    }

    AtomSet::from_unsorted(hits)
}

fn cell_of(
    p: [f32; 3],
    inv: f32,
    pbc: bool,
    box_diag: [f32; 3],
    ncells: [i32; 3],
) -> (i32, i32, i32) {
    let axis = |k: usize| {
        if pbc {
            // Fold the coordinate into the box first, so an atom that drifted outside still lands
            // in the cell its periodic image occupies.
            let l = box_diag[k];
            let wrapped = p[k] - l * (p[k] / l).floor();
            ((wrapped * inv).floor() as i32).rem_euclid(ncells[k])
        } else {
            (p[k] * inv).floor() as i32
        }
    };
    (axis(0), axis(1), axis(2))
}

fn dist2(p: [f32; 3], q: [f32; 3], box_diag: [f32; 3], pbc: bool) -> f32 {
    let mut sum = 0.0;
    for k in 0..3 {
        let mut d = q[k] - p[k];
        if pbc {
            // Minimum image, rectangular box.
            let l = box_diag[k];
            d -= l * (d / l).round();
        }
        sum += d * d;
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line_of_carbons(xs: &[f32], box_x: f32) -> Structure {
        Structure {
            names: vec!["C".into(); xs.len()],
            resnames: vec!["X".into(); xs.len()],
            resids: vec![1; xs.len()],
            elements: vec!["C".into(); xs.len()],
            positions: xs.iter().map(|&x| [x, 0.0, 0.0]).collect(),
            box_diag: [box_x, box_x, box_x],
            triclinic: false,
        }
    }

    #[test]
    fn finds_atoms_inside_the_radius() {
        // Atoms at x = 0, 0.3, 0.6, 2.0
        let s = line_of_carbons(&[0.0, 0.3, 0.6, 2.0], 10.0);
        let of = AtomSet::from_unsorted(vec![1]);
        assert_eq!(within(&s, &of, 0.5, false).as_slice(), [1, 2]);
        assert_eq!(within(&s, &of, 0.65, false).as_slice(), [1, 2, 3]);
    }

    #[test]
    fn the_selection_is_inside_its_own_neighbourhood() {
        let s = line_of_carbons(&[0.0, 5.0], 10.0);
        let of = AtomSet::from_unsorted(vec![1]);
        // Nothing else is near, but atom 1 is at distance 0 from itself.
        assert_eq!(within(&s, &of, 0.1, false).as_slice(), [1]);
    }

    #[test]
    fn an_empty_selection_has_an_empty_neighbourhood() {
        let s = line_of_carbons(&[0.0, 1.0], 10.0);
        assert!(within(&s, &AtomSet::empty(), 5.0, false).is_empty());
    }

    /// Across the periodic boundary the two atoms are 0.2 nm apart, not 9.8.
    #[test]
    fn pbc_wraps_across_the_box() {
        let s = line_of_carbons(&[0.1, 9.9], 10.0);
        let of = AtomSet::from_unsorted(vec![1]);
        assert_eq!(
            within(&s, &of, 0.5, false).as_slice(),
            [1],
            "without --pbc they are far apart"
        );
        assert_eq!(
            within(&s, &of, 0.5, true).as_slice(),
            [1, 2],
            "with --pbc they are neighbours"
        );
    }

    #[test]
    fn a_zero_radius_selects_only_coincident_atoms() {
        let s = line_of_carbons(&[0.0, 0.0, 1.0], 10.0);
        let of = AtomSet::from_unsorted(vec![1]);
        assert_eq!(within(&s, &of, 0.0, false).as_slice(), [1, 2]);
    }

    #[test]
    fn a_large_radius_selects_everything() {
        let s = line_of_carbons(&[0.0, 1.0, 2.0, 3.0], 10.0);
        let of = AtomSet::from_unsorted(vec![1]);
        assert_eq!(within(&s, &of, 100.0, false).as_slice(), [1, 2, 3, 4]);
    }
}
