pub mod preprocess;

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crate::error::{NdxError, Result};
use crate::model::AtomId;

use preprocess::Preprocessor;

/// One `[ moleculetype ]`: the template a `[ molecules ]` line instantiates.
#[derive(Clone, Debug, Default)]
struct Template {
    name: String,
    natoms: usize,
    types: Vec<String>,
    charges: Vec<f32>,
    masses: Vec<f32>,
    /// 1-based within the molecule.
    bonds: Vec<(u32, u32)>,
}

/// A whole system's topology, expanded across `[ molecules ]` so every index is a global,
/// 1-based atom number — the same numbering the `.ndx` uses.
#[derive(Clone, Debug, Default)]
pub struct Topology {
    pub natoms: u32,

    // Interned: a 40k-atom system has a handful of distinct types and molecule names.
    type_ids: Vec<u32>,
    pub type_names: Vec<String>,
    mol_ids: Vec<u32>,
    pub mol_names: Vec<String>,
    /// 0-based instance number: the 3rd copy of SOL is 2. This is what `splitmol` splits on.
    pub mol_instance: Vec<u32>,

    pub charges: Vec<f32>,
    pub masses: Vec<f32>,

    /// Global 1-based bond pairs, `a < b`, deduplicated.
    pub bonds: Vec<(AtomId, AtomId)>,

    /// `#include`s that could not be found. Harmless if they only held parameters.
    pub missing_includes: Vec<String>,
}

impl Topology {
    fn idx(&self, a: AtomId) -> Option<usize> {
        let i = (a as usize).checked_sub(1)?;
        (i < self.type_ids.len()).then_some(i)
    }

    pub fn atomtype(&self, a: AtomId) -> Option<&str> {
        self.idx(a)
            .map(|i| self.type_names[self.type_ids[i] as usize].as_str())
    }

    pub fn molname(&self, a: AtomId) -> Option<&str> {
        self.idx(a)
            .map(|i| self.mol_names[self.mol_ids[i] as usize].as_str())
    }

    pub fn instance(&self, a: AtomId) -> Option<u32> {
        self.idx(a).map(|i| self.mol_instance[i])
    }

    pub fn charge(&self, a: AtomId) -> Option<f32> {
        self.idx(a).map(|i| self.charges[i])
    }

    pub fn mass(&self, a: AtomId) -> Option<f32> {
        self.idx(a).map(|i| self.masses[i])
    }

    pub fn parse_path(
        path: &Path,
        defines: &[String],
        include_dirs: &[PathBuf],
    ) -> Result<Self> {
        let mut pre = Preprocessor::new(defines, include_dirs);
        let expanded = pre.expand_path(path)?;
        let mut top = Self::parse_expanded(&expanded, &path.display().to_string())?;
        top.missing_includes = std::mem::take(&mut pre.missing);
        Ok(top)
    }

    /// Parse an already-preprocessed topology (no `#` directives left).
    pub fn parse_expanded(content: &str, origin: &str) -> Result<Self> {
        let (templates, instances, inter_bonds) = parse_sections(content);

        // Expand: walk the [ molecules ] list, laying each copy down after the last.
        let mut top = Topology::default();
        let mut intern_type = Interner::default();
        let mut intern_mol = Interner::default();
        let mut offset: u32 = 0;
        // Instance numbers run per species, not per `[ molecules ]` line: a topology may list
        // `SOL 100` and later `SOL 50`, and those are copies 0..149 of the same molecule.
        let mut next_instance: HashMap<u32, u32> = HashMap::new();

        for (mol_name, copies) in &instances {
            let Some(t) = templates.iter().find(|t| &t.name == mol_name) else {
                return Err(NdxError::UnknownMoleculeType {
                    name: mol_name.clone(),
                    known: templates.iter().map(|t| t.name.clone()).collect(),
                });
            };
            let mol_id = intern_mol.intern(&t.name);
            let type_ids: Vec<u32> = t.types.iter().map(|s| intern_type.intern(s)).collect();

            for _ in 0..*copies {
                let instance = next_instance.entry(mol_id).or_insert(0);
                for ((&ty, &charge), &mass) in type_ids
                    .iter()
                    .zip(&t.charges)
                    .zip(&t.masses)
                {
                    top.type_ids.push(ty);
                    top.charges.push(charge);
                    top.masses.push(mass);
                    top.mol_ids.push(mol_id);
                    top.mol_instance.push(*instance);
                }
                for &(ai, aj) in &t.bonds {
                    // Template bonds are 1-based within the molecule.
                    top.bonds.push((ai + offset, aj + offset));
                }
                *instance += 1;
                offset += t.natoms as u32;
            }
        }

        // `[ intermolecular_interactions ]` bonds are already global 1-based numbers.
        top.bonds.extend(inter_bonds);

        top.natoms = offset;
        top.type_names = intern_type.names;
        top.mol_names = intern_mol.names;

        // Normalize: a < b, sorted, unique. Both the graph builder and the tests want this.
        for b in &mut top.bonds {
            if b.0 > b.1 {
                std::mem::swap(&mut b.0, &mut b.1);
            }
        }
        top.bonds.sort_unstable();
        top.bonds.dedup();

        if let Some(&(_, hi)) = top.bonds.iter().max_by_key(|b| b.1)
            && hi > top.natoms
        {
            return Err(NdxError::Other(format!(
                "{origin}: a bond references atom {hi}, but the topology expands to only {} atoms",
                top.natoms
            )));
        }

        Ok(top)
    }
}

#[derive(Default)]
struct Interner {
    names: Vec<String>,
}

impl Interner {
    fn intern(&mut self, s: &str) -> u32 {
        if let Some(i) = self.names.iter().position(|n| n == s) {
            return i as u32;
        }
        self.names.push(s.to_string());
        (self.names.len() - 1) as u32
    }
}

type Instances = Vec<(String, usize)>;

fn parse_sections(content: &str) -> (Vec<Template>, Instances, Vec<(AtomId, AtomId)>) {
    let mut templates: Vec<Template> = Vec::new();
    let mut instances: Instances = Vec::new();
    let mut inter_bonds: Vec<(AtomId, AtomId)> = Vec::new();

    let mut section = String::new();
    let mut intermolecular = false;
    let mut current: Option<Template> = None;
    // A `[ moleculetype ]` header is followed by one data line holding its name.
    let mut want_name = false;

    for raw in content.lines() {
        let line = strip_comment(raw);
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        if let Some(name) = section_header(line) {
            if name == "moleculetype" {
                if let Some(t) = current.take() {
                    templates.push(t);
                }
                current = Some(Template::default());
                want_name = true;
            }
            if name == "intermolecular_interactions" {
                // These bonds use global atom numbers, so they must not go into a template.
                if let Some(t) = current.take() {
                    templates.push(t);
                }
                intermolecular = true;
            }
            section = name;
            continue;
        }

        if want_name && section == "moleculetype" {
            if let Some(t) = current.as_mut() {
                t.name = line.split_whitespace().next().unwrap_or("").to_string();
            }
            want_name = false;
            continue;
        }

        match section.as_str() {
            "molecules" => {
                let f: Vec<&str> = line.split_whitespace().collect();
                if f.len() >= 2
                    && let Ok(n) = f[1].parse::<usize>()
                {
                    instances.push((f[0].to_string(), n));
                }
            }

            "atoms" => {
                if let Some(t) = current.as_mut()
                    && let Some(a) = parse_atom(line)
                {
                    t.types.push(a.0);
                    t.charges.push(a.1);
                    t.masses.push(a.2);
                    t.natoms += 1;
                }
            }

            // Bonds proper, plus the two sections that are *also* bonds as far as connectivity
            // goes. Without them a rigid water has no bonds at all, and a system run with
            // `constraints = h-bonds` loses every X-H.
            "bonds" | "constraints" => {
                let f: Vec<&str> = line.split_whitespace().collect();
                if f.len() >= 2
                    && let (Ok(ai), Ok(aj)) = (f[0].parse::<u32>(), f[1].parse::<u32>())
                    && ai > 0
                    && aj > 0
                {
                    if intermolecular {
                        inter_bonds.push((ai, aj));
                    } else if let Some(t) = current.as_mut() {
                        t.bonds.push((ai, aj));
                    }
                }
            }

            // `[ settles ]` is `OW funct doh dhh`: the named atom plus the next two are one rigid
            // water. GROMACS never writes those two bonds anywhere else.
            "settles" => {
                if let Some(t) = current.as_mut() {
                    let f: Vec<&str> = line.split_whitespace().collect();
                    if !f.is_empty()
                        && let Ok(ow) = f[0].parse::<u32>()
                        && ow > 0
                    {
                        t.bonds.push((ow, ow + 1));
                        t.bonds.push((ow, ow + 2));
                    }
                }
            }

            _ => {}
        }
    }

    if let Some(t) = current {
        templates.push(t);
    }

    (templates, instances, inter_bonds)
}

/// `nr type resnr residue atom cgnr [charge] [mass]`
///
/// Charge and mass are optional — they default to the atom type's values, which live in a force
/// field we may not have. Returns them as 0.0 / NaN-free defaults; `Structure` only uses the mass
/// when it actually identifies an element.
fn parse_atom(line: &str) -> Option<(String, f32, f32)> {
    let f: Vec<&str> = line.split_whitespace().collect();
    if f.len() < 6 {
        return None;
    }
    f[0].parse::<u32>().ok()?; // the atom number; position is what actually orders them
    let atom_type = f[1].to_string();
    let charge = f.get(6).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let mass = f.get(7).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    Some((atom_type, charge, mass))
}

fn section_header(line: &str) -> Option<String> {
    let inner = line.strip_prefix('[')?.strip_suffix(']')?;
    Some(inner.trim().to_ascii_lowercase())
}

fn strip_comment(line: &str) -> &str {
    match line.find(';') {
        Some(i) => &line[..i],
        None => line,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOP: &str = "\
[ moleculetype ]
ETH   3

[ atoms ]
    1  opls_135  1   ETH    C1     1   -0.180   12.011
    2  opls_140  1   ETH    H11    1    0.060    1.008
    3  opls_140  1   ETH    H12    1    0.060    1.008

[ bonds ]
    1     2     1
    1     3     1

[ moleculetype ]
SOL   2

[ atoms ]
    1  opls_116  1   SOL    OW     1   -0.834   15.999
    2  opls_117  1   SOL    HW1    1    0.417    1.008
    3  opls_117  1   SOL    HW2    1    0.417    1.008

[ settles ]
    1     1   0.09572  0.15139

[ molecules ]
ETH   2
SOL   3
";

    fn top() -> Topology {
        Topology::parse_expanded(TOP, "t").unwrap()
    }

    #[test]
    fn molecules_are_expanded_with_the_right_offsets() {
        let t = top();
        assert_eq!(t.natoms, 2 * 3 + 3 * 3);
        assert_eq!(t.mol_names, ["ETH", "SOL"]);

        // The second ETH copy starts at atom 4.
        assert_eq!(t.atomtype(1), Some("opls_135"));
        assert_eq!(t.atomtype(4), Some("opls_135"));
        assert_eq!(t.molname(4), Some("ETH"));
        // ... and the waters follow it.
        assert_eq!(t.molname(7), Some("SOL"));
        assert_eq!(t.atomtype(7), Some("opls_116"));
    }

    /// Bonds are template-local; each copy has to be shifted by the atoms laid down before it.
    #[test]
    fn bonds_are_offset_per_copy() {
        let t = top();
        assert!(t.bonds.contains(&(1, 2)), "first ETH");
        assert!(t.bonds.contains(&(4, 5)), "second ETH, shifted by 3");
        assert!(!t.bonds.contains(&(3, 4)), "copies must not be joined together");
    }

    /// Rigid water has no [ bonds ] at all — its connectivity is implied by [ settles ]. Missing
    /// this means every water in the system comes out with zero bonds.
    #[test]
    fn settles_become_bonds() {
        let t = top();
        // First water is atoms 7,8,9: OW-HW1 and OW-HW2.
        assert!(t.bonds.contains(&(7, 8)));
        assert!(t.bonds.contains(&(7, 9)));
        assert!(!t.bonds.contains(&(8, 9)), "the two H are not bonded to each other");
        // 2 bonds per ETH x 2, plus 2 per water x 3.
        assert_eq!(t.bonds.len(), 2 * 2 + 2 * 3);
    }

    #[test]
    fn instances_are_numbered_per_species() {
        let t = top();
        assert_eq!(t.instance(1), Some(0), "first ETH");
        assert_eq!(t.instance(4), Some(1), "second ETH");
        assert_eq!(t.instance(7), Some(0), "first SOL restarts at 0");
        assert_eq!(t.instance(10), Some(1));
        assert_eq!(t.instance(13), Some(2));
    }

    /// A topology may name the same species in more than one [ molecules ] block.
    #[test]
    fn a_species_split_across_blocks_keeps_counting() {
        let src = TOP.replace("SOL   3", "SOL   2\nSOL   1");
        let t = Topology::parse_expanded(&src, "t").unwrap();
        assert_eq!(t.natoms, 15);
        assert_eq!(t.instance(7), Some(0));
        assert_eq!(t.instance(13), Some(2), "the third water is copy 2, not copy 0");
    }

    #[test]
    fn constraints_count_as_bonds() {
        let src = "\
[ moleculetype ]
M  3
[ atoms ]
 1 t 1 M A 1 0.0 1.0
 2 t 1 M B 1 0.0 1.0
[ constraints ]
 1 2 1 0.1
[ molecules ]
M 1
";
        let t = Topology::parse_expanded(src, "t").unwrap();
        assert_eq!(t.bonds, [(1, 2)]);
    }

    #[test]
    fn charges_and_masses_are_read() {
        let t = top();
        assert!((t.mass(1).unwrap() - 12.011).abs() < 1e-3);
        assert!((t.charge(1).unwrap() + 0.180).abs() < 1e-3);
    }

    /// The atom line's charge and mass are optional (they default to the atom type's).
    #[test]
    fn an_atom_line_without_charge_or_mass_still_parses() {
        let src = "\
[ moleculetype ]
M  3
[ atoms ]
 1 opls_1 1 M A 1
[ molecules ]
M 1
";
        let t = Topology::parse_expanded(src, "t").unwrap();
        assert_eq!(t.natoms, 1);
        assert_eq!(t.atomtype(1), Some("opls_1"));
        assert_eq!(t.mass(1), Some(0.0));
    }

    /// The definition is usually in an #include we could not find, so say so.
    #[test]
    fn an_unknown_molecule_type_is_a_clear_error() {
        let src = "[ molecules ]\nSOL 100\n";
        match Topology::parse_expanded(src, "t") {
            Err(NdxError::UnknownMoleculeType { name, .. }) => assert_eq!(name, "SOL"),
            other => panic!("expected UnknownMoleculeType, got {other:?}"),
        }
        let msg = Topology::parse_expanded(src, "t").unwrap_err().to_string();
        assert!(msg.contains("#include"), "{msg}");
        assert!(msg.contains("GMXLIB"), "{msg}");
    }

    /// These use global atom numbers, not per-molecule ones.
    #[test]
    fn intermolecular_bonds_are_global() {
        let src = format!("{TOP}\n[ intermolecular_interactions ]\n[ bonds ]\n 1 7 1\n");
        let t = Topology::parse_expanded(&src, "t").unwrap();
        assert!(t.bonds.contains(&(1, 7)), "an ETH atom bonded to a water atom");
    }

    #[test]
    fn comments_are_stripped() {
        let src = "\
[ moleculetype ] ; the molecule
M  3
[ atoms ]
 1 t 1 M A 1 0.0 1.0 ; a comment
[ molecules ]
M 1  ; one copy
";
        assert_eq!(Topology::parse_expanded(src, "t").unwrap().natoms, 1);
    }
}
