//! Tests against the real system in `test_data/`: 512,328 atoms, a bundled fiber in solvent.
//!
//! That directory is ~35 MB, so it is not something every checkout should be required to have.
//! Each test **skips** when it is absent rather than failing, which keeps `cargo test` green on a
//! bare clone while still exercising the real thing when the data is there.
//!
//! What makes this system worth testing against, beyond its size:
//!   - the `.gro` carries velocities, so its atom lines run past the coordinate columns
//!   - `topo.top` `#include`s `MCH.itp`, and `mixed_hbond.itp` behind `#ifdef INTER`
//!   - `mixed_hbond.itp` is an `[ intermolecular_interactions ]` block: global atom numbers
//!   - the `.ndx` covers only the fiber (atoms 1..40920), not the whole box

mod common;

use std::path::{Path, PathBuf};

use common::*;
use predicates::prelude::*;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("test_data")
        .join(name)
}

/// `None` when the data is not checked out, which every test treats as "skip".
fn system() -> Option<(String, String, String)> {
    let (ndx, gro, top) = (data("index.ndx"), data("input.gro"), data("topo.top"));
    if !(ndx.is_file() && gro.is_file() && top.is_file()) {
        eprintln!("skipping: test_data/ is not present");
        return None;
    }
    Some((
        ndx.to_str().unwrap().into(),
        gro.to_str().unwrap().into(),
        top.to_str().unwrap().into(),
    ))
}

macro_rules! real {
    ($ndx:ident, $gro:ident, $top:ident) => {
        let Some(($ndx, $gro, $top)) = system() else {
            return;
        };
    };
}

fn select(ndx: &str, expr: &str, flags: &[&str]) -> Vec<u32> {
    let mut args = vec!["select".to_string(), ndx.to_string(), expr.to_string()];
    args.extend(flags.iter().map(|s| (*s).to_string()));
    args.push("--keep-only".into());
    args.push("-q".into());

    let out = ndxed()
        .args(&args)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    parse_out(&out).groups[0].atoms.clone()
}

/// The structure and the topology have to agree, or every element and type lookup is off by
/// however far they drifted apart. They do here — 512,328 atoms both ways.
#[test]
fn the_structure_and_topology_describe_the_same_system() {
    real!(ndx, gro, top);

    ndxed()
        .args(["info", &ndx, "-s", &gro, "-p", &top, "-q"])
        .assert()
        .success()
        .stdout(predicate::str::contains("structure: 512328 atoms"))
        .stdout(predicate::str::contains("topology:  512328 atoms"))
        .stdout(predicate::str::contains("BAE, BAS, MCH"))
        .stdout(predicate::str::contains("(from the topology)"))
        // No "do not match" warning: the two agree.
        .stderr(predicate::str::contains("do not match").not());
}

/// A `.gro` written with velocities has six numbers after the atom name, not three. Reading it by
/// fixed columns is what keeps that from corrupting the coordinates.
#[test]
fn a_gro_with_velocities_parses() {
    real!(ndx, gro, _top);
    // If the velocity columns had been mistaken for coordinates, the atoms would be scattered and
    // this shell would come out the wrong size.
    let near = select(&ndx, "within 0.2 of atomid 1", &["-s", &gro]);
    assert!(
        (2..40).contains(&near.len()),
        "0.2 nm around one atom should hold a handful of neighbours, got {}",
        near.len()
    );
    assert!(near.contains(&1), "an atom is inside its own shell");
}

/// The query this was all for.
#[test]
fn the_hydrogens_bonded_to_the_fiber() {
    real!(ndx, gro, top);
    let h = select(&ndx, "element H & bonded Fiber1", &["-s", &gro, "-p", &top]);
    assert_eq!(h.len(), 18864);

    // Every one of them must be in Fiber1 and be a hydrogen — check by construction.
    let fiber = select(&ndx, "Fiber1", &["-s", &gro, "-p", &top]);
    let all_h = select(&ndx, "element H & Fiber1", &["-s", &gro, "-p", &top]);
    assert!(h.iter().all(|a| fiber.contains(a)));
    assert!(h.iter().all(|a| all_h.contains(a)));
}

/// `#ifdef INTER` guards the intermolecular hydrogen bonds. Honouring it is not optional: with the
/// conditional ignored, 264 bonds that the simulation does not have would appear in every answer.
#[test]
fn ifdef_gates_the_intermolecular_bonds() {
    real!(ndx, gro, top);

    let count = |flags: &[&str]| {
        let mut args = vec!["info", "-p", &top, "-q"];
        args.extend_from_slice(flags);
        let out = ndxed().args(&args).assert().success().get_output().stdout.clone();
        let s = String::from_utf8(out).unwrap();
        s.lines()
            .find(|l| l.starts_with("bonds:"))
            .and_then(|l| l.split_whitespace().nth(1))
            .and_then(|n| n.parse::<usize>().ok())
            .expect("a bond count")
    };

    let without = count(&[]);
    let with = count(&["-D", "INTER"]);
    // mixed_hbond.itp holds exactly 264 [ bonds ] lines.
    assert_eq!(with - without, 264, "-D INTER must add exactly the h-bonds");

    // And they must actually be wired into the graph, at their global atom numbers: the first
    // pair in the file is 116 - 15781.
    let plain = select(&ndx, "bonded atomid 116", &["-s", &gro, "-p", &top]);
    let inter = select(
        &ndx,
        "bonded atomid 116",
        &["-s", &gro, "-p", &top, "-D", "INTER"],
    );
    assert!(!plain.contains(&15781));
    assert!(
        inter.contains(&15781),
        "the [ intermolecular_interactions ] bond 116-15781 should appear with -D INTER"
    );
    assert_eq!(inter.len(), plain.len() + 1);
}

/// With no topology, bonds fall back to a distance estimate. It should find essentially the same
/// covalent bonds — that is the whole basis for the fallback being useful.
#[test]
fn distance_inference_broadly_agrees_with_the_topology_on_a_real_system() {
    real!(ndx, gro, top);

    let from_top = select(&ndx, "element H & bonded Core1", &["-s", &gro, "-p", &top]);
    let from_dist = select(&ndx, "element H & bonded Core1", &["-s", &gro]);

    let agree = from_top.iter().filter(|a| from_dist.contains(a)).count();
    let union = from_top.len() + from_dist.len() - agree;
    let jaccard = agree as f64 / union as f64;
    assert!(
        jaccard > 0.95,
        "distance-inferred bonds should closely match the topology's: \
         {agree}/{union} = {jaccard:.3} (top {}, distance {})",
        from_top.len(),
        from_dist.len()
    );
}

#[test]
fn topology_predicates_on_a_real_system() {
    real!(ndx, gro, top);
    let p = ["-p", &top];
    let _ = &gro;

    // The .ndx only covers the fiber; the molecule predicate sees the whole 512k-atom system, so
    // intersecting with a group is how you scope it.
    let mch = select(&ndx, "molecule MCH & Fiber1", &p);
    assert!(mch.is_empty(), "the solvent is not part of Fiber1");

    let fiber_mols = select(&ndx, "(molecule BAE | molecule BAS) & Fiber1", &p);
    assert_eq!(fiber_mols.len(), 33264, "Fiber1 is entirely fiber molecules");
}

/// The .ndx stops at atom 40920 but the system has 512,328. That is legal — an index file need not
/// cover everything — and `!` must complement against the *system*, not the index file.
#[test]
fn the_universe_comes_from_the_structure_not_the_index_file() {
    real!(ndx, gro, top);

    // Without a structure, `!Fiber1` can only complement against the groups it can see.
    let guessed = select(&ndx, "!Fiber1", &[]);
    // With one, it complements against all 512,328 atoms.
    let exact = select(&ndx, "!Fiber1", &["-s", &gro, "-p", &top]);

    assert_eq!(exact.len(), 512328 - 33264);
    assert!(exact.len() > guessed.len());
}

#[test]
fn the_real_index_file_round_trips_byte_for_byte() {
    real!(ndx, _gro, _top);
    let out = ndxed()
        .args(["fmt", &ndx, "-q"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(out, bytes(Path::new(&ndx)));
}
