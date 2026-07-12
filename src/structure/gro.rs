//! GROMACS `.gro` coordinate files.
//!
//! Fixed-column format, and it must be read as such — a residue name and an atom name can run
//! together with no space between them, so whitespace splitting silently corrupts such files.
//!
//! ```text
//! MD of 2 waters, t= 0.0
//!     6
//!     1WATER  OW1    1   0.126   1.624   1.679
//! |----|----|----|----|-------|-------|-------|
//!  0..5 5..10 10..15 15..20 20..28  28..36  36..44
//! ```
//!
//! Line 1 is a title, line 2 the atom count, then one line per atom, then the box vectors.
//! Coordinates are in **nm**.

use std::path::Path;

use crate::error::{NdxError, Result};

#[derive(Clone, Debug, PartialEq)]
pub struct GroAtom {
    pub resid: i32,
    pub resname: String,
    pub name: String,
    /// The `.gro`'s own atom number. Not trusted for indexing — it wraps at 99999.
    pub serial: i64,
    /// nm
    pub pos: [f32; 3],
}

#[derive(Clone, Debug, Default)]
pub struct GroFile {
    pub title: String,
    pub atoms: Vec<GroAtom>,
    /// Box vectors in nm: the three diagonal elements, then the off-diagonals if triclinic.
    pub box_diag: [f32; 3],
    pub triclinic: bool,
}

impl GroFile {
    pub fn parse(content: &str, origin: &str) -> Result<Self> {
        let mut lines = content.lines();

        let title = lines.next().unwrap_or_default().trim_end().to_string();

        let count_line = lines
            .next()
            .ok_or_else(|| bad(origin, 2, "expected an atom count on line 2"))?;
        let declared: usize = count_line.trim().parse().map_err(|_| {
            bad(
                origin,
                2,
                &format!("expected an atom count, got {:?}", count_line.trim()),
            )
        })?;

        let mut atoms = Vec::with_capacity(declared);
        let mut box_diag = [0.0f32; 3];
        let mut triclinic = false;

        for (i, line) in lines.enumerate() {
            let lineno = i + 3;

            // The atom count is authoritative: whatever follows the declared atoms is the box.
            if atoms.len() == declared {
                if line.trim().is_empty() {
                    continue;
                }
                let vals: Vec<f32> = line
                    .split_ascii_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if vals.len() < 3 {
                    return Err(bad(
                        origin,
                        lineno,
                        "expected 3 (or 9) box vector components",
                    ));
                }
                box_diag = [vals[0], vals[1], vals[2]];
                // A triclinic box carries 6 more components; a non-zero one means the simple
                // rectangular minimum-image convention no longer applies.
                triclinic = vals.len() > 3 && vals[3..].iter().any(|v| v.abs() > 1e-6);
                break;
            }

            let atom = parse_atom_line(line).ok_or_else(|| {
                bad(
                    origin,
                    lineno,
                    &format!("malformed atom line: {:?}", line.trim_end()),
                )
            })?;
            atoms.push(atom);
        }

        if atoms.len() != declared {
            return Err(NdxError::Other(format!(
                "{origin}: the header declares {declared} atoms but the file has {}",
                atoms.len()
            )));
        }

        Ok(GroFile {
            title,
            atoms,
            box_diag,
            triclinic,
        })
    }

    pub fn parse_path(p: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(p).map_err(|e| NdxError::io(p, e))?;
        Self::parse(&content, &p.display().to_string())
    }
}

/// Slice by byte column, tolerating a line that stops early.
fn col(line: &str, from: usize, to: usize) -> Option<&str> {
    let b = line.as_bytes();
    if b.len() < from {
        return None;
    }
    let to = to.min(b.len());
    // The format is ASCII; guard anyway so a stray multi-byte char cannot panic.
    line.get(from..to)
}

fn parse_atom_line(line: &str) -> Option<GroAtom> {
    let resid: i32 = col(line, 0, 5)?.trim().parse().ok()?;
    let resname = col(line, 5, 10)?.trim().to_string();
    let name = col(line, 10, 15)?.trim().to_string();
    // The serial wraps to 0 after 99999 in big systems, so it is read but never used to index.
    let serial: i64 = col(line, 15, 20)?.trim().parse().ok()?;
    let x: f32 = col(line, 20, 28)?.trim().parse().ok()?;
    let y: f32 = col(line, 28, 36)?.trim().parse().ok()?;
    let z: f32 = col(line, 36, 44)?.trim().parse().ok()?;

    if name.is_empty() {
        return None;
    }

    Some(GroAtom {
        resid,
        resname,
        name,
        serial,
        pos: [x, y, z],
    })
}

fn bad(origin: &str, line: usize, msg: &str) -> NdxError {
    NdxError::Other(format!("{origin}:{line}: {msg}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
MD of 2 waters
    6
    1SOL     OW    1   0.230   0.628   0.113
    1SOL    HW1    2   0.137   0.626   0.150
    1SOL    HW2    3   0.231   0.589   0.021
    2SOL     OW    4   1.230   1.628   1.113
    2SOL    HW1    5   1.137   1.626   1.150
    2SOL    HW2    6   1.231   1.589   1.021
   1.82060   1.82060   1.82060
";

    #[test]
    fn parses_a_gro() {
        let g = GroFile::parse(SAMPLE, "t").unwrap();
        assert_eq!(g.title, "MD of 2 waters");
        assert_eq!(g.atoms.len(), 6);
        assert_eq!(g.atoms[0].resname, "SOL");
        assert_eq!(g.atoms[0].name, "OW");
        assert_eq!(g.atoms[0].resid, 1);
        assert_eq!(g.atoms[3].resid, 2);
        assert!((g.atoms[0].pos[0] - 0.230).abs() < 1e-5);
        assert!((g.box_diag[0] - 1.8206).abs() < 1e-4);
        assert!(!g.triclinic);
    }

    /// The whole reason for reading fixed columns: with a 5-digit residue number and 5-character
    /// names, these fields touch with no space anywhere. Whitespace splitting reads one token here.
    #[test]
    fn handles_fields_that_run_together() {
        let src = "t\n1\n99999SOLUTEATOM99999   1.000   2.000   3.000\n1 1 1\n";
        let g = GroFile::parse(src, "t").unwrap();
        assert_eq!(g.atoms[0].resid, 99999);
        assert_eq!(g.atoms[0].resname, "SOLUT");
        assert_eq!(g.atoms[0].name, "EATOM");
        assert_eq!(g.atoms[0].serial, 99999);
        assert_eq!(g.atoms[0].pos, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn detects_a_triclinic_box() {
        let src = "t\n1\n    1SOL     OW    1   0.100   0.200   0.300\n2.0 2.0 2.0 0.0 0.0 1.0 0.0 0.0 0.0\n";
        let g = GroFile::parse(src, "t").unwrap();
        assert!(g.triclinic);
    }

    #[test]
    fn the_declared_count_is_authoritative() {
        // Two atoms declared, so the third line is the box, not an atom.
        let src = "t\n2\n    1A       X    1   0.100   0.200   0.300\n    1A       Y    2   0.100   0.200   0.300\n3.0 3.0 3.0\n";
        let g = GroFile::parse(src, "t").unwrap();
        assert_eq!(g.atoms.len(), 2);
        assert_eq!(g.box_diag, [3.0, 3.0, 3.0]);
    }

    #[test]
    fn a_short_file_is_an_error() {
        let src = "t\n5\n    1A       X    1   0.100   0.200   0.300\n";
        let err = GroFile::parse(src, "t").unwrap_err().to_string();
        assert!(err.contains("declares 5 atoms"), "{err}");
    }

    #[test]
    fn a_bad_count_line_is_an_error() {
        assert!(GroFile::parse("t\nnope\n", "t").is_err());
    }
}
