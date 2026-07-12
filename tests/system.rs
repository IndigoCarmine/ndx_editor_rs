//! Structure- and topology-aware selection, end to end.
//!
//! `system.gro` / `system.top` / `system.ndx` describe one small system:
//!   atoms  1..18  two ethanol molecules (ETH), 9 atoms each: C1 H11 H12 H13 C2 H21 H22 O HO
//!   atoms 19..27  three rigid waters (SOL), bonded only via `[ settles ]`
//!   atom     28   a sodium ion (NA), its own residue

mod common;

use common::*;
use predicates::prelude::*;

fn ndx() -> String {
    fixture("system.ndx").to_str().unwrap().to_string()
}
fn gro() -> String {
    fixture("system.gro").to_str().unwrap().to_string()
}
fn top() -> String {
    fixture("system.top").to_str().unwrap().to_string()
}

/// Run `select` and return the atoms of the resulting group.
fn select(expr: &str, flags: &[&str]) -> Vec<u32> {
    let mut args = vec!["select".to_string(), ndx(), expr.to_string()];
    args.extend(flags.iter().map(|s| (*s).to_string()));
    args.push("--keep-only".into());

    let out = ndxed()
        .args(&args)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    parse_out(&out).groups[0].atoms.clone()
}

// ------------------------------------------------------------ the target query

/// The whole point of reading a `.gro` and a `.top`: the hydrogens attached to a group.
#[test]
fn the_hydrogens_bonded_to_a_group() {
    let h = select("element H & bonded ETH", &["-s", &gro(), "-p", &top()]);
    // Six hydrogens per ethanol (H11 H12 H13 H21 H22 HO), two molecules.
    assert_eq!(h, [2, 3, 4, 6, 7, 9, 11, 12, 13, 15, 16, 18]);
}

/// The same answer from coordinates alone, which is what makes `-s` without `-p` useful.
#[test]
fn distance_inferred_bonds_agree_with_the_topology() {
    let from_top = select("element H & bonded ETH", &["-s", &gro(), "-p", &top()]);
    let from_dist = select("element H & bonded ETH", &["-s", &gro()]);
    assert_eq!(from_dist, from_top);
}

#[test]
fn estimating_bonds_says_that_it_is_estimating() {
    ndxed()
        .args(["select", &ndx(), "bonded ETH", "-s", &gro()])
        .assert()
        .success()
        .stderr(predicate::str::contains("estimated from interatomic distances"))
        .stderr(predicate::str::contains("-p topol.top"));
}

// ---------------------------------------------------------------------- bonds

#[test]
fn bonded_is_the_neighbourhood_of_the_selection() {
    let flags = ["-s", &gro(), "-p", &top()];

    // Ethanol is a closed molecule: everything bonded to it is already in it.
    assert_eq!(select("bonded ETH", &flags), (1..=18).collect::<Vec<u32>>());
    assert!(select("bonded ETH & !ETH", &flags).is_empty());

    // Two bonds out from the ethanol oxygens reaches C2, its hydrogens, C1, and the hydroxyl H.
    assert_eq!(
        select("bonded 2 of name O", &flags),
        [1, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18]
    );
}

/// Rigid water carries no `[ bonds ]` — only `[ settles ]`. Missing that would leave every water
/// in the system with no connectivity at all.
#[test]
fn settles_gives_water_its_bonds() {
    let flags = ["-s", &gro(), "-p", &top()];
    // OW of the first water is atom 19; its hydrogens are 20 and 21.
    assert_eq!(select("bonded name OW", &flags), [20, 21, 23, 24, 26, 27]);
}

#[test]
fn info_reports_the_bond_count_and_its_source() {
    ndxed()
        .args(["info", &ndx(), "-s", &gro(), "-p", &top()])
        .assert()
        .success()
        // 8 bonds per ethanol x2, plus 2 per water x3, from [ settles ].
        .stdout(predicate::str::contains("bonds:     22 (from the topology)"))
        .stdout(predicate::str::contains("28 atoms"))
        .stdout(predicate::str::contains("ETH, SOL, NA"));
}

// ----------------------------------------------------------------- predicates

#[test]
fn structure_predicates() {
    let s = ["-s", &gro()];
    assert_eq!(select("resname SOL", &s), (19..=27).collect::<Vec<u32>>());
    assert_eq!(select("resid 1-2", &s), (1..=18).collect::<Vec<u32>>());
    assert_eq!(select("name OW", &s), [19, 22, 25]);
    assert_eq!(select("name HW1 HW2", &s), [20, 21, 23, 24, 26, 27]);
    assert_eq!(select("name H*", &s).len(), 18, "6 per ethanol + 2 per water");
    assert_eq!(select("element C", &s), [1, 5, 10, 14]);
    assert_eq!(select("element O", &s), [8, 17, 19, 22, 25]);
}

#[test]
fn topology_predicates() {
    let p = ["-p", &top()];
    assert_eq!(select("molecule SOL", &p), (19..=27).collect::<Vec<u32>>());
    assert_eq!(select("molecule ETH", &p), (1..=18).collect::<Vec<u32>>());
    // opls_140 is the aliphatic hydrogen type: 5 per ethanol (the hydroxyl H is opls_155).
    assert_eq!(select("type opls_140", &p), [2, 3, 4, 6, 7, 11, 12, 13, 15, 16]);
}

/// `CA` in a protein is a carbon, but a lone `NA` residue really is sodium. The mass settles it.
#[test]
fn a_monatomic_ion_gets_its_real_element() {
    assert_eq!(select("element NA", &["-s", &gro(), "-p", &top()]), [28]);
}

#[test]
fn predicates_compose_with_group_algebra_and_wildcards() {
    let flags = ["-s", &gro(), "-p", &top()];
    // Water hydrogens only: `element H` intersected with a group from the .ndx.
    assert_eq!(select("element H & SOL", &flags), [20, 21, 23, 24, 26, 27]);
    // Everything that is not water or ion.
    assert_eq!(select("System & !(SOL | NA)", &flags), (1..=18).collect::<Vec<u32>>());
}

// -------------------------------------------------------------------- within

#[test]
fn within_finds_nearby_atoms() {
    let s = ["-s", &gro()];

    // Each water's hydrogens sit ~0.1 nm from its oxygen, so a 0.2 nm shell around the three OW
    // atoms is exactly the three waters — and nothing else, since they are 0.4 nm apart.
    assert_eq!(select("within 0.2 of name OW", &s), (19..=27).collect::<Vec<u32>>());

    // The ion is isolated: the nearest water is 1.58 nm away.
    assert_eq!(select("within 0.4 of NA", &s), [28]);
    assert_eq!(select("within 1.0 of NA", &s), [28]);
    assert!(select("within 2.0 of NA", &s).len() > 1, "a wide enough shell reaches the waters");

    // A selection is inside its own neighbourhood; `& !X` is the strictly-outside form.
    assert_eq!(select("within 0.2 of name OW & !name OW", &s), [20, 21, 23, 24, 26, 27]);
}

// --------------------------------------------------------------------- split

#[test]
fn split_by_molecule() {
    let out = ndxed()
        .args([
            "split", &ndx(), "ETH", "--by", "molecule", "-p", &top(), "--replace",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let f = parse_out(&out);
    let n = names(&f);
    assert!(n.contains(&"ETH_ETH1".to_string()), "{n:?}");
    assert!(n.contains(&"ETH_ETH2".to_string()), "{n:?}");
    let eth1 = f.groups.iter().find(|g| g.name == "ETH_ETH1").unwrap();
    assert_eq!(eth1.atoms, (1..=9).collect::<Vec<u32>>());
}

#[test]
fn split_by_residue() {
    let out = ndxed()
        .args(["split", &ndx(), "SOL", "--by", "residue", "-s", &gro()])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let f = parse_out(&out);
    let n = names(&f);
    // The waters are residues 3, 4 and 5.
    assert!(n.contains(&"SOL_SOL3".to_string()), "{n:?}");
    assert!(n.contains(&"SOL_SOL5".to_string()), "{n:?}");
    let sol3 = f.groups.iter().find(|g| g.name == "SOL_SOL3").unwrap();
    assert_eq!(sol3.atoms, [19, 20, 21]);
}

#[test]
fn split_by_residue_without_a_structure_says_so() {
    ndxed()
        .args(["split", &ndx(), "SOL", "--by", "residue"])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("-s conf.gro"));
}

// ------------------------------------------------- missing files and mismatches

/// A force-field `#include` we cannot resolve holds parameters, not molecules — carry on, but say
/// so, because a *missing molecule* would otherwise look like a mystery.
#[test]
fn a_missing_include_warns_but_still_works() {
    ndxed()
        .args(["info", "-p", &top()])
        .assert()
        .success()
        .stderr(predicate::str::contains("oplsaa.ff/forcefield.itp"))
        .stdout(predicate::str::contains("28 atoms"));
}

#[test]
fn a_structure_that_does_not_match_the_index_is_caught() {
    let dir = tempfile::tempdir().unwrap();
    let big = dir.path().join("big.ndx");
    std::fs::write(&big, "[ Big ]\n 9999\n").unwrap();

    ndxed()
        .args(["info", big.to_str().unwrap(), "-s", &gro()])
        .assert()
        .success()
        .stderr(predicate::str::contains("references atom 9999"))
        .stderr(predicate::str::contains("do not match"));
}

#[test]
fn a_structure_and_topology_of_different_sizes_are_caught() {
    let dir = tempfile::tempdir().unwrap();
    let small = dir.path().join("small.gro");
    std::fs::write(
        &small,
        "one atom\n    1\n    1SOL     OW    1   0.100   0.200   0.300\n1 1 1\n",
    )
    .unwrap();

    ndxed()
        .args(["info", "-s", small.to_str().unwrap(), "-p", &top()])
        .assert()
        .success()
        .stderr(predicate::str::contains("1 atoms but the topology expands to 28"));
}

/// Without `-s`/`-p` the expressions still parse, and the error names the flag that would work.
#[test]
fn without_a_system_the_error_names_the_flag() {
    ndxed()
        .args(["select", &ndx(), "element H & bonded ETH"])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("`element` needs a structure file"))
        .stderr(predicate::str::contains("-s conf.gro"));

    ndxed()
        .args(["select", &ndx(), "type opls_140"])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("`type` needs a topology"));

    // With only a topology there are no coordinates, so `within` still cannot run.
    ndxed()
        .args(["select", &ndx(), "within 0.5 of 1", "-p", &top()])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("needs atom coordinates"));
}

/// A .gro has no chain column; point at the thing GROMACS actually has.
#[test]
fn chain_points_at_molecule_instead() {
    ndxed()
        .args(["select", &ndx(), "chain A", "-s", &gro()])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("no chain column"))
        .stderr(predicate::str::contains("molecule <name>"));
}

// ----------------------------------------------------------- the interactive editor

#[test]
fn the_editor_loads_the_system_and_uses_it() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.ndx");

    let assert = ndxed()
        .args([
            "edit", &ndx(), "-s", &gro(), "-p", &top(), "-o",
            out.to_str().unwrap(),
        ])
        .write_stdin("element H & bonded ETH\nname 4 H_on_ETH\nq\n")
        .assert()
        .success();

    let stdout = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert!(stdout.contains("bonds:     22"), "{stdout}");
    assert!(stdout.contains("12 atoms"), "{stdout}");

    let f = ndx_editor::parse_path(&out).unwrap();
    let g = f.groups.iter().find(|g| g.name == "H_on_ETH").unwrap();
    assert_eq!(g.atoms, [2, 3, 4, 6, 7, 9, 11, 12, 13, 15, 16, 18]);
}
