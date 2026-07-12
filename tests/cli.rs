mod common;

use common::*;
use predicates::prelude::*;

// ---------------------------------------------------------------- round-trip

/// The load-bearing invariant: a gmx-formatted file survives a read/write cycle byte for byte.
#[test]
fn gmx_format_round_trips_byte_for_byte() {
    for name in ["small.ndx", "spaces.ndx", "dup_names.ndx", "no_system.ndx"] {
        let path = fixture(name);
        let out = ndxed()
            .args(["fmt", path.to_str().unwrap()])
            .assert()
            .success()
            .get_output()
            .stdout
            .clone();
        assert_eq!(out, bytes(&path), "{name} did not round-trip");
    }
}

#[test]
fn a_messy_file_parses_and_normalizes() {
    let out = ndxed()
        .args(["fmt", fixture("messy.ndx").to_str().unwrap()])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let f = parse_out(&out);
    assert_eq!(names(&f), ["System", "Empty", "Unpadded"]);
    assert_eq!(f.groups[0].atoms, [1, 2, 3, 4, 5]);
    assert!(f.groups[1].atoms.is_empty());
}

// ---------------------------------------------------------------------- list

#[test]
fn list_shows_groups_and_counts() {
    ndxed()
        .args(["list", fixture("small.ndx").to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains("System"))
        .stdout(predicate::str::contains("50 atoms"))
        .stdout(predicate::str::contains("Protein"));
}

#[test]
fn list_long_flags_unsorted_groups() {
    ndxed()
        .args(["list", "--long", fixture("unsorted.ndx").to_str().unwrap()])
        .assert()
        .success()
        .stdout(predicate::str::contains("unsorted"));
}

#[test]
fn list_json_is_parseable_shape() {
    ndxed()
        .args([
            "list",
            "--format",
            "json",
            fixture("small.ndx").to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains(r#""name": "Protein""#))
        .stdout(predicate::str::contains(r#""atoms": 20"#));
}

#[test]
fn list_reads_stdin() {
    ndxed()
        .args(["list", "-"])
        .write_stdin(String::from_utf8(bytes(&fixture("small.ndx"))).unwrap())
        .assert()
        .success()
        .stdout(predicate::str::contains("System"));
}

// -------------------------------------------------------------------- select

#[test]
fn select_make_ndx_expression() {
    let out = ndxed()
        .args(["select", fixture("small.ndx").to_str().unwrap(), "1 & !2"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let f = parse_out(&out);
    assert_eq!(f.len(), 4);
    assert_eq!(f.groups[3].name, "Protein_&_!SOL");
    assert_eq!(f.groups[3].atoms, (1..=20).collect::<Vec<u32>>());
}

#[test]
fn select_by_name_with_spaces() {
    let out = ndxed()
        .args([
            "select",
            fixture("spaces.ndx").to_str().unwrap(),
            "C-alpha | Water and ions",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let f = parse_out(&out);
    assert_eq!(f.groups[3].atoms, [2, 5, 8, 9, 10]);
}

#[test]
fn select_can_be_piped_into_list() {
    let selected = ndxed()
        .args(["select", fixture("small.ndx").to_str().unwrap(), "1 | 2"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    ndxed()
        .args(["list", "-"])
        .write_stdin(selected)
        .assert()
        .success()
        .stdout(predicate::str::contains("Protein_|_SOL"));
}

#[test]
fn select_keep_only_drops_the_rest() {
    let out = ndxed()
        .args([
            "select",
            fixture("small.ndx").to_str().unwrap(),
            "1 & 2",
            "--keep-only",
            "--name",
            "Overlap",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let f = parse_out(&out);
    assert_eq!(names(&f), ["Overlap"]);
}

/// The rewrite in `expr::eval`: `A & !B` must not need a universe, so no warning here even though
/// this file has no `System` group.
#[test]
fn and_not_does_not_warn_about_the_universe() {
    ndxed()
        .args(["select", fixture("no_system.ndx").to_str().unwrap(), "0 & !1"])
        .assert()
        .success()
        .stderr(predicate::str::contains("universe").not());
}

/// A bare `!` genuinely does need one, and must say so.
#[test]
fn a_bare_complement_warns_when_the_universe_is_guessed() {
    ndxed()
        .args(["select", fixture("no_system.ndx").to_str().unwrap(), "!1"])
        .assert()
        .success()
        .stderr(predicate::str::contains("union of all groups"))
        .stderr(predicate::str::contains("--natoms"));
}

#[test]
fn natoms_makes_the_complement_exact_and_silent() {
    let out = ndxed()
        .args([
            "select",
            fixture("no_system.ndx").to_str().unwrap(),
            "!0",
            "--natoms",
            "8",
            "--name",
            "NotProtein",
        ])
        .assert()
        .success()
        .stderr(predicate::str::contains("warning").not())
        .get_output()
        .stdout
        .clone();
    let f = parse_out(&out);
    assert_eq!(f.groups[2].atoms, [4, 5, 6, 7, 8]);
}

#[test]
fn empty_selection_warns_but_still_creates_the_group() {
    let out = ndxed()
        .args(["select", fixture("small.ndx").to_str().unwrap(), "1 & !1"])
        .assert()
        .success()
        .stderr(predicate::str::contains("no atoms"))
        .get_output()
        .stdout
        .clone();
    assert_eq!(parse_out(&out).len(), 4);
}

#[test]
fn error_on_empty_opts_out() {
    ndxed()
        .args([
            "select",
            fixture("small.ndx").to_str().unwrap(),
            "1 & !1",
            "--error-on-empty",
        ])
        .assert()
        .code(3);
}

#[test]
fn duplicate_names_are_ambiguous_unless_numbered() {
    let f = fixture("dup_names.ndx");
    ndxed()
        .args(["select", f.to_str().unwrap(), "SOL"])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("ambiguous"))
        .stderr(predicate::str::contains("SOL#1"));

    let out = ndxed()
        .args(["select", f.to_str().unwrap(), "SOL#2"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(parse_out(&out).groups[3].atoms, [3, 4]);
}

#[test]
fn an_unknown_group_suggests_a_near_miss() {
    ndxed()
        .args(["select", fixture("small.ndx").to_str().unwrap(), "Protien"])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("did you mean"))
        .stderr(predicate::str::contains("Protein"));
}

#[test]
fn a_syntax_error_points_a_caret_at_the_problem() {
    ndxed()
        .args(["select", fixture("small.ndx").to_str().unwrap(), "0 & "])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("missing an operand"))
        .stderr(predicate::str::contains("^"));
}

// ------------------------------------- structure-dependent features, with no system loaded

/// These need a `.gro` / `.top`. Without one the failure must name the flag that would have
/// supplied it, rather than reading as a syntax error. (See tests/system.rs for them working.)
#[test]
fn structure_expressions_report_the_flag_they_need() {
    let f = fixture("small.ndx");

    ndxed()
        .args(["select", f.to_str().unwrap(), "element H & bonded Protein"])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("`element` needs a structure file"))
        .stderr(predicate::str::contains("-s conf.gro"))
        .stderr(predicate::str::contains("^"));

    ndxed()
        .args(["select", f.to_str().unwrap(), "bonded Protein"])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("`bonded` needs bond information"))
        .stderr(predicate::str::contains("topol.top"));

    ndxed()
        .args(["select", f.to_str().unwrap(), "within 0.5 of Protein"])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("`within` needs atom coordinates"));

    ndxed()
        .args(["select", f.to_str().unwrap(), "type OW"])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("`type` needs a topology"));
}

// -------------------------------------------------- wildcards in group names

/// `Fiber*` stands for the union of every group it matches.
#[test]
fn a_glob_selects_the_union_of_its_matches() {
    let out = ndxed()
        .args(["select", fixture("globby.ndx").to_str().unwrap(), "Fiber*"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let f = parse_out(&out);
    assert_eq!(f.groups[5].name, "Fiber*");
    assert_eq!(f.groups[5].atoms, [1, 2, 3, 4]);
}

#[test]
fn globs_compose_with_the_operators() {
    let f = fixture("globby.ndx");
    let sel = |expr: &str| {
        let out = ndxed()
            .args(["select", f.to_str().unwrap(), expr, "--keep-only"])
            .assert()
            .success()
            .get_output()
            .stdout
            .clone();
        parse_out(&out).groups[0].atoms.clone()
    };
    // Fiber1={1,2} Fiber2={3,4} fiberA={5} AlkylA={6} OA={7}
    assert_eq!(sel("Fiber* | fiber*"), [1, 2, 3, 4, 5]);
    assert_eq!(sel("Fiber* & !Fiber2"), [1, 2]);
    assert_eq!(sel("*A"), [5, 6, 7]);
    assert_eq!(sel("O?"), [7]);
}

#[test]
fn del_and_keep_take_globs() {
    let f = fixture("globby.ndx");

    let out = ndxed()
        .args(["del", f.to_str().unwrap(), "Fiber*"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(names(&parse_out(&out)), ["fiberA", "AlkylA", "OA"]);

    let out = ndxed()
        .args(["keep", f.to_str().unwrap(), "*A"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(names(&parse_out(&out)), ["fiberA", "AlkylA", "OA"]);
}

/// A pattern that matches nothing must be an error, not a silent no-op that quietly writes an
/// unchanged file.
#[test]
fn a_glob_matching_nothing_is_an_error() {
    ndxed()
        .args(["del", fixture("globby.ndx").to_str().unwrap(), "Nope*"])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("no group name matches"));
}

/// rename/split/--universe act on one group, so a glob that matches several must be refused
/// rather than picking one arbitrarily.
#[test]
fn commands_that_need_one_group_refuse_an_ambiguous_glob() {
    ndxed()
        .args([
            "rename",
            fixture("globby.ndx").to_str().unwrap(),
            "Fiber*",
            "X",
        ])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("matches 2 groups"))
        .stderr(predicate::str::contains("Fiber1, Fiber2"));

    // A glob with exactly one match is fine.
    let out = ndxed()
        .args([
            "rename",
            fixture("globby.ndx").to_str().unwrap(),
            "Alkyl*",
            "Chains",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert!(names(&parse_out(&out)).contains(&"Chains".to_string()));
}

/// Quoting turns the glob off, so a group genuinely named `O*` is still reachable.
#[test]
fn quoting_makes_a_wildcard_literal() {
    let dir = tempfile::tempdir().unwrap();
    let star = dir.path().join("star.ndx");
    std::fs::write(&star, "[ O* ]\n   1    2\n[ OA ]\n   3\n").unwrap();

    let sel = |expr: &str| {
        let out = ndxed()
            .args(["select", star.to_str().unwrap(), expr, "--keep-only"])
            .assert()
            .success()
            .get_output()
            .stdout
            .clone();
        parse_out(&out).groups[0].atoms.clone()
    };
    assert_eq!(sel("O*"), [1, 2, 3], "bare: a pattern, matching both groups");
    assert_eq!(sel("\"O*\""), [1, 2], "quoted: the group actually named O*");
}

// ------------------------------------------------------------ simple mutators

#[test]
fn rename() {
    let out = ndxed()
        .args([
            "rename",
            fixture("small.ndx").to_str().unwrap(),
            "Protein",
            "Prot",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(names(&parse_out(&out)), ["System", "Prot", "SOL"]);
}

#[test]
fn del_and_keep() {
    let f = fixture("small.ndx");

    let out = ndxed()
        .args(["del", f.to_str().unwrap(), "0"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(names(&parse_out(&out)), ["Protein", "SOL"]);

    // `keep` also reorders, since it takes the groups in the order given.
    let out = ndxed()
        .args(["keep", f.to_str().unwrap(), "SOL", "0"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(names(&parse_out(&out)), ["SOL", "System"]);
}

#[test]
fn head_and_tail_take_atoms_in_file_order() {
    let f = fixture("unsorted.ndx");

    let out = ndxed()
        .args(["head", f.to_str().unwrap(), "Backwards", "2"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    // The group is 5 4 3 2 1, so its *first* two atoms are 5 and 4 — not 1 and 2.
    assert_eq!(parse_out(&out).groups[2].atoms, [5, 4]);

    let out = ndxed()
        .args(["head", f.to_str().unwrap(), "Backwards", "2", "--tail"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(parse_out(&out).groups[2].atoms, [2, 1]);
}

#[test]
fn head_strict_refuses_to_clamp() {
    ndxed()
        .args([
            "head",
            fixture("small.ndx").to_str().unwrap(),
            "Protein",
            "999",
            "--strict",
        ])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("only 20 atom"));
}

#[test]
fn atoms_from_literal_ranges() {
    let out = ndxed()
        .args([
            "atoms",
            fixture("small.ndx").to_str().unwrap(),
            "1-3,7",
            "--name",
            "Mine",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let f = parse_out(&out);
    assert_eq!(f.groups[3].name, "Mine");
    assert_eq!(f.groups[3].atoms, [1, 2, 3, 7]);
}

// --------------------------------------------------------------------- split

#[test]
fn split_parts_size_and_at() {
    let f = fixture("small.ndx");

    let out = ndxed()
        .args(["split", f.to_str().unwrap(), "Protein", "--parts", "3"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let g = parse_out(&out);
    // 20 atoms into 3 -> 7, 7, 6
    assert_eq!(names(&g)[3..], ["Protein_1", "Protein_2", "Protein_3"]);
    assert_eq!(
        [g.groups[3].len(), g.groups[4].len(), g.groups[5].len()],
        [7, 7, 6]
    );

    let out = ndxed()
        .args(["split", f.to_str().unwrap(), "Protein", "--size", "8"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let g = parse_out(&out);
    assert_eq!(
        [g.groups[3].len(), g.groups[4].len(), g.groups[5].len()],
        [8, 8, 4]
    );

    let out = ndxed()
        .args([
            "split",
            f.to_str().unwrap(),
            "Protein",
            "--at",
            "5,15",
            "--replace",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let g = parse_out(&out);
    assert_eq!(names(&g), ["System", "SOL", "Protein_1", "Protein_2", "Protein_3"]);
    assert_eq!(g.groups[2].atoms, [1, 2, 3, 4, 5]);
    assert_eq!(g.groups[4].atoms, [16, 17, 18, 19, 20]);
}

#[test]
fn split_preserves_atom_order() {
    let out = ndxed()
        .args([
            "split",
            fixture("unsorted.ndx").to_str().unwrap(),
            "Backwards",
            "--size",
            "2",
            "--replace",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let f = parse_out(&out);
    assert_eq!(f.groups[1].atoms, [5, 4]);
    assert_eq!(f.groups[2].atoms, [3, 2]);
    assert_eq!(f.groups[3].atoms, [1]);
}

#[test]
fn split_needs_exactly_one_mode() {
    ndxed()
        .args([
            "split",
            fixture("small.ndx").to_str().unwrap(),
            "Protein",
            "--parts",
            "2",
            "--size",
            "2",
        ])
        .assert()
        .code(2);
}

// --------------------------------------------------------------------- merge

#[test]
fn merge_concatenates_with_each_conflict_policy() {
    let (a, b) = (fixture("small.ndx"), fixture("other.ndx"));
    let run = |policy: &str| {
        let out = ndxed()
            .args([
                "merge",
                a.to_str().unwrap(),
                b.to_str().unwrap(),
                "--on-conflict",
                policy,
            ])
            .assert()
            .success()
            .get_output()
            .stdout
            .clone();
        parse_out(&out)
    };

    let f = run("keep-both");
    assert_eq!(names(&f), ["System", "Protein", "SOL", "Protein", "Ligand"]);
    assert_eq!(f.groups[3].atoms, [100, 101]);

    let f = run("rename");
    assert_eq!(
        names(&f),
        ["System", "Protein", "SOL", "Protein_2", "Ligand"]
    );

    let f = run("skip");
    assert_eq!(names(&f), ["System", "Protein", "SOL", "Ligand"]);
    assert_eq!(f.groups[1].atoms.len(), 20, "the first Protein survives");

    let f = run("replace");
    assert_eq!(names(&f), ["System", "Protein", "SOL", "Ligand"]);
    assert_eq!(f.groups[1].atoms, [100, 101], "the last Protein wins");
}

#[test]
fn merge_can_prefix_by_file() {
    let out = ndxed()
        .args([
            "merge",
            fixture("small.ndx").to_str().unwrap(),
            fixture("other.ndx").to_str().unwrap(),
            "--prefix-file",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(
        names(&parse_out(&out)),
        [
            "small:System",
            "small:Protein",
            "small:SOL",
            "other:Protein",
            "other:Ligand"
        ]
    );
}

// ---------------------------------------------------------------------- diff

#[test]
fn diff_of_a_file_with_itself_exits_zero() {
    let f = fixture("small.ndx");
    ndxed()
        .args(["diff", f.to_str().unwrap(), f.to_str().unwrap()])
        .assert()
        .code(0)
        .stdout(predicate::str::contains("no differences"));
}

#[test]
fn diff_exits_one_when_files_differ() {
    ndxed()
        .args([
            "diff",
            fixture("small.ndx").to_str().unwrap(),
            fixture("other.ndx").to_str().unwrap(),
        ])
        .assert()
        .code(1)
        .stdout(predicate::str::contains("- System"))
        .stdout(predicate::str::contains("+ Ligand"))
        .stdout(predicate::str::contains("~ Protein"));
}

#[test]
fn diff_no_exit_code_stays_quiet_for_scripts() {
    ndxed()
        .args([
            "diff",
            fixture("small.ndx").to_str().unwrap(),
            fixture("other.ndx").to_str().unwrap(),
            "--no-exit-code",
        ])
        .assert()
        .code(0);
}

#[test]
fn diff_atoms_itemizes_the_change() {
    // no_system.ndx has Protein = 1 2 3. Change it to 2 3 99: one atom gone, one gained.
    let dir = tempfile::tempdir().unwrap();
    let edited = dir.path().join("edited.ndx");
    std::fs::write(&edited, "[ Protein ]\n   2    3   99\n[ SOL ]\n   4    5\n").unwrap();

    ndxed()
        .args([
            "diff",
            fixture("no_system.ndx").to_str().unwrap(),
            edited.to_str().unwrap(),
            "--atoms",
        ])
        .assert()
        .code(1)
        .stdout(predicate::str::contains("~ Protein (3 -> 3 atoms; +1 -1)"))
        // The itemized lines: atom 1 was removed, atom 99 added.
        .stdout(predicate::str::contains("- 1"))
        .stdout(predicate::str::contains("+ 99"));
}

#[test]
fn diff_detects_a_pure_reorder() {
    let dir = tempfile::tempdir().unwrap();
    let sorted = dir.path().join("sorted.ndx");
    ndxed()
        .args([
            "fmt",
            fixture("unsorted.ndx").to_str().unwrap(),
            "--sort",
            "-o",
            sorted.to_str().unwrap(),
        ])
        .assert()
        .success();

    ndxed()
        .args([
            "diff",
            fixture("unsorted.ndx").to_str().unwrap(),
            sorted.to_str().unwrap(),
        ])
        .assert()
        .code(1)
        .stdout(predicate::str::contains("same atoms, different order"));
}

// ------------------------------------------------------- errors and file safety

/// Nothing may ever modify the file it was given.
#[test]
fn no_command_touches_its_input() {
    let f = fixture("small.ndx");
    let before = bytes(&f);
    let p = f.to_str().unwrap();

    for args in [
        vec!["select", p, "1 & !2"],
        vec!["rename", p, "0", "Renamed"],
        vec!["del", p, "0"],
        vec!["keep", p, "1"],
        vec!["head", p, "1", "3"],
        vec!["atoms", p, "1-5"],
        vec!["split", p, "1", "--parts", "2"],
        vec!["fmt", p, "--sort", "--dedup"],
    ] {
        ndxed().args(&args).assert().success();
        assert_eq!(bytes(&f), before, "{args:?} modified its input");
    }
}

#[test]
fn writing_over_the_input_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let target = dir.path().join("work.ndx");
    std::fs::write(&target, bytes(&fixture("small.ndx"))).unwrap();
    let before = bytes(&target);

    ndxed()
        .args([
            "select",
            target.to_str().unwrap(),
            "1 & !2",
            "-o",
            target.to_str().unwrap(),
        ])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("refusing to overwrite"));

    assert_eq!(bytes(&target), before, "the input must be untouched");
}

#[test]
fn force_overwrite_opts_in() {
    let dir = tempfile::tempdir().unwrap();
    let target = dir.path().join("work.ndx");
    std::fs::write(&target, bytes(&fixture("small.ndx"))).unwrap();

    ndxed()
        .args([
            "select",
            target.to_str().unwrap(),
            "1 & !2",
            "-o",
            target.to_str().unwrap(),
            "--force-overwrite",
        ])
        .assert()
        .success();

    let f = ndx_editor::parse_path(&target).unwrap();
    assert_eq!(f.len(), 4);
}

#[test]
fn output_to_a_file_writes_it() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.ndx");
    ndxed()
        .args([
            "select",
            fixture("small.ndx").to_str().unwrap(),
            "1 | 2",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::is_empty())
        .stderr(predicate::str::contains("wrote 4 group(s)"));

    assert_eq!(ndx_editor::parse_path(&out).unwrap().len(), 4);
}

#[test]
fn a_malformed_file_is_rejected_with_a_location() {
    ndxed()
        .args(["list", fixture("malformed_zero.ndx").to_str().unwrap()])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("1-based"))
        .stderr(predicate::str::contains(":2:"));

    ndxed()
        .args(["list", fixture("malformed_early.ndx").to_str().unwrap()])
        .assert()
        .code(3)
        .stderr(predicate::str::contains("before any"));
}

#[test]
fn a_missing_file_is_an_io_error() {
    ndxed()
        .args(["list", "/nonexistent/nope.ndx"])
        .assert()
        .code(4);
}

#[test]
fn bad_flags_are_a_usage_error() {
    ndxed().args(["list", "--bogus"]).assert().code(2);
}

// -------------------------------------------------------------- completions

#[test]
fn completions_generate_for_every_shell() {
    for shell in ["bash", "zsh", "fish", "elvish", "powershell"] {
        ndxed()
            .args(["completions", shell])
            .assert()
            .success()
            .stdout(predicate::str::contains("ndxed"));
    }
}
