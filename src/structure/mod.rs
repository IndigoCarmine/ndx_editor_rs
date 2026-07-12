pub mod gro;

use std::path::Path;

use crate::error::Result;
use crate::model::AtomId;

pub use gro::GroFile;

/// Per-atom identity and coordinates, indexed by 1-based [`AtomId`] to match the `.ndx`.
///
/// Struct-of-arrays: a predicate sweeps one field over every atom, so keeping the fields in
/// separate contiguous vectors is both simpler and faster than a vector of structs.
#[derive(Clone, Debug, Default)]
pub struct Structure {
    pub names: Vec<String>,
    pub resnames: Vec<String>,
    pub resids: Vec<i32>,
    /// Element symbol: `H`, `C`, `CL`, ... Uppercase. `X` when it could not be worked out.
    pub elements: Vec<String>,
    /// nm
    pub positions: Vec<[f32; 3]>,
    /// Box diagonal in nm, for `--pbc`.
    pub box_diag: [f32; 3],
    pub triclinic: bool,
}

impl Structure {
    pub fn natoms(&self) -> u32 {
        self.names.len() as u32
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// `.ndx` atom ids are 1-based; the arrays are 0-based.
    fn idx(&self, a: AtomId) -> Option<usize> {
        let i = (a as usize).checked_sub(1)?;
        (i < self.names.len()).then_some(i)
    }

    pub fn name(&self, a: AtomId) -> Option<&str> {
        self.idx(a).map(|i| self.names[i].as_str())
    }

    pub fn resname(&self, a: AtomId) -> Option<&str> {
        self.idx(a).map(|i| self.resnames[i].as_str())
    }

    pub fn resid(&self, a: AtomId) -> Option<i32> {
        self.idx(a).map(|i| self.resids[i])
    }

    pub fn element(&self, a: AtomId) -> Option<&str> {
        self.idx(a).map(|i| self.elements[i].as_str())
    }

    pub fn pos(&self, a: AtomId) -> Option<[f32; 3]> {
        self.idx(a).map(|i| self.positions[i])
    }

    pub fn from_gro(gro: &GroFile) -> Self {
        let n = gro.atoms.len();
        let mut s = Structure {
            names: Vec::with_capacity(n),
            resnames: Vec::with_capacity(n),
            resids: Vec::with_capacity(n),
            elements: Vec::with_capacity(n),
            positions: Vec::with_capacity(n),
            box_diag: gro.box_diag,
            triclinic: gro.triclinic,
        };
        for a in &gro.atoms {
            s.elements.push(element_from_name(&a.name, &a.resname));
            s.names.push(a.name.clone());
            s.resnames.push(a.resname.clone());
            s.resids.push(a.resid);
            s.positions.push(a.pos);
        }
        s
    }

    pub fn from_gro_path(p: &Path) -> Result<Self> {
        Ok(Self::from_gro(&GroFile::parse_path(p)?))
    }

    /// Refine the elements using masses from a topology.
    ///
    /// A mass identifies an element far more reliably than a name does — the name `CA` is an alpha
    /// carbon in a protein but a calcium ion on its own, and a mass of 12.011 settles it. Masses
    /// that match nothing (coarse-grained beads; hydrogen-mass repartitioning, which moves H to
    /// ~3 amu) are left to the name-based guess rather than forced into a wrong element.
    pub fn refine_elements_from_masses(&mut self, masses: &[f32]) -> usize {
        let mut refined = 0;
        for (i, &m) in masses.iter().enumerate().take(self.elements.len()) {
            if let Some(sym) = element_from_mass(m)
                && self.elements[i] != sym
            {
                self.elements[i] = sym.to_string();
                refined += 1;
            }
        }
        refined
    }
}

/// Elements we can name, with their standard atomic masses.
const MASSES: &[(&str, f32)] = &[
    ("H", 1.008),
    ("C", 12.011),
    ("N", 14.007),
    ("O", 15.999),
    ("F", 18.998),
    ("NA", 22.990),
    ("MG", 24.305),
    ("P", 30.974),
    ("S", 32.060),
    ("CL", 35.450),
    ("K", 39.098),
    ("CA", 40.078),
    ("FE", 55.845),
    ("CU", 63.546),
    ("ZN", 65.380),
    ("BR", 79.904),
    ("I", 126.904),
];

/// Two-letter elements that also occur as ordinary atom names inside a residue (`CA` is both an
/// alpha carbon and calcium), so they only count as the element when the atom stands alone as its
/// own residue — the way monatomic ions are written.
const AMBIGUOUS_TWO_LETTER: &[&str] = &["CA", "NA", "CL", "CU", "FE", "MG", "ZN", "K", "BR"];

fn element_from_mass(m: f32) -> Option<&'static str> {
    MASSES
        .iter()
        .find(|(_, mass)| (m - mass).abs() < 0.15)
        .map(|(sym, _)| *sym)
}

/// Guess an element from a GROMACS atom name.
///
/// The rule that matters: **the first alphabetic character is the element**, so `CA`, `CB`, `CG1`
/// are all carbon and `HW1`, `1HB` are all hydrogen. That is right for the overwhelming majority
/// of atoms in a force field. The exception is a monatomic ion, where the atom is its own residue
/// (`resname CL`, `atom CL`) — there the whole name is the element.
pub fn element_from_name(name: &str, resname: &str) -> String {
    let name = name.trim();
    // PDB-style names may lead with a digit: `1HB` is a hydrogen.
    let alpha: String = name
        .chars()
        .skip_while(|c| !c.is_ascii_alphabetic())
        .take_while(|c| c.is_ascii_alphabetic())
        .collect::<String>()
        .to_ascii_uppercase();

    if alpha.is_empty() {
        return "X".to_string();
    }

    // A monatomic ion: the atom name is the residue name, so the two-letter reading is right.
    let res = resname.trim().to_ascii_uppercase();
    if alpha == res && AMBIGUOUS_TWO_LETTER.contains(&alpha.as_str()) {
        return alpha;
    }
    if MASSES.iter().any(|(s, _)| *s == alpha) && alpha.len() == 2 && alpha == res {
        return alpha;
    }

    alpha.chars().next().expect("alpha is non-empty").to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_first_letter_is_the_element() {
        assert_eq!(element_from_name("CA", "ALA"), "C", "alpha carbon");
        assert_eq!(element_from_name("CB", "ALA"), "C");
        assert_eq!(element_from_name("CG1", "VAL"), "C");
        assert_eq!(element_from_name("OW", "SOL"), "O");
        assert_eq!(element_from_name("HW1", "SOL"), "H");
        assert_eq!(element_from_name("N", "ALA"), "N");
        assert_eq!(element_from_name("SD", "MET"), "S");
    }

    #[test]
    fn pdb_style_leading_digits() {
        assert_eq!(element_from_name("1HB", "ALA"), "H");
        assert_eq!(element_from_name("2HG1", "ILE"), "H");
    }

    /// The case a naive first-letter rule gets wrong.
    #[test]
    fn a_monatomic_ion_keeps_both_letters() {
        assert_eq!(element_from_name("CA", "CA"), "CA", "calcium ion");
        assert_eq!(element_from_name("CL", "CL"), "CL");
        assert_eq!(element_from_name("NA", "NA"), "NA");
        assert_eq!(element_from_name("ZN", "ZN"), "ZN");
        // ... but the same name inside a real residue is still a carbon.
        assert_eq!(element_from_name("CA", "LYS"), "C");
    }

    #[test]
    fn unknown_names_fall_back() {
        assert_eq!(element_from_name("", "X"), "X");
        assert_eq!(element_from_name("123", "X"), "X");
    }

    #[test]
    fn masses_identify_elements() {
        assert_eq!(element_from_mass(1.008), Some("H"));
        assert_eq!(element_from_mass(12.011), Some("C"));
        assert_eq!(element_from_mass(12.01), Some("C"), "force fields round");
        assert_eq!(element_from_mass(15.9994), Some("O"));
        assert_eq!(element_from_mass(40.08), Some("CA"));
        // Coarse-grained bead: no element, and we must not invent one.
        assert_eq!(element_from_mass(72.0), None);
        // Hydrogen-mass repartitioning moves H to ~3 amu; better to fall back to the name than
        // to claim it is helium.
        assert_eq!(element_from_mass(3.024), None);
    }

    #[test]
    fn masses_override_an_ambiguous_name() {
        let mut s = Structure {
            names: vec!["CA".into()],
            resnames: vec!["ALA".into()],
            resids: vec![1],
            elements: vec![element_from_name("CA", "ALA")],
            positions: vec![[0.0; 3]],
            ..Default::default()
        };
        assert_eq!(s.elements[0], "C");
        // A topology saying 40.078 means it really was calcium.
        assert_eq!(s.refine_elements_from_masses(&[40.078]), 1);
        assert_eq!(s.elements[0], "CA");
    }

    #[test]
    fn an_unrecognized_mass_leaves_the_name_guess_alone() {
        let mut s = Structure {
            names: vec!["BB".into()],
            resnames: vec!["MARTINI".into()],
            resids: vec![1],
            elements: vec!["B".into()],
            positions: vec![[0.0; 3]],
            ..Default::default()
        };
        assert_eq!(s.refine_elements_from_masses(&[72.0]), 0);
        assert_eq!(s.elements[0], "B");
    }

    #[test]
    fn structure_is_indexed_from_one() {
        let gro = GroFile::parse(
            "t\n2\n    1SOL     OW    1   0.100   0.200   0.300\n    1SOL    HW1    2   0.400   0.500   0.600\n2 2 2\n",
            "t",
        )
        .unwrap();
        let s = Structure::from_gro(&gro);
        assert_eq!(s.natoms(), 2);
        assert_eq!(s.name(1), Some("OW"));
        assert_eq!(s.element(1), Some("O"));
        assert_eq!(s.name(2), Some("HW1"));
        assert_eq!(s.element(2), Some("H"));
        assert_eq!(s.name(0), None, "there is no atom 0");
        assert_eq!(s.name(3), None);
    }
}
