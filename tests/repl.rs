//! The REPL, driven through stdin.
//!
//! The non-tty path is a plain line loop, so these exercise the real command dispatch, the real
//! session, and the real writer — with no pty anywhere.

mod common;

use common::*;
use predicates::prelude::*;

/// Feed a script to `ndxed edit` and return (stdout, exit code).
fn script(input: &str, script: &str) -> (String, i32) {
    let out = ndxed()
        .args(["edit", fixture(input).to_str().unwrap()])
        .write_stdin(script.to_string())
        .assert()
        .get_output()
        .clone();
    (
        String::from_utf8(out.stdout).unwrap(),
        out.status.code().unwrap(),
    )
}

#[test]
fn it_greets_with_the_group_table() {
    let (out, code) = script("small.ndx", "q!\n");
    assert_eq!(code, 0);
    assert!(out.contains("Reading index file"));
    assert!(out.contains("System"));
    assert!(out.contains("50 atoms"));
    assert!(out.contains("Type `help`"));
}

#[test]
fn a_bare_expression_adds_a_group_and_q_saves_it() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.ndx");

    let (out, code) = script(
        "small.ndx",
        &format!("l\n0 & !1\nq {}\n", out_path.display()),
    );
    assert_eq!(code, 0);
    assert!(out.contains("System_&_!Protein"), "{out}");
    assert!(out.contains("wrote 4 group(s)"));

    // Assert on the file that was actually written.
    let f = ndx_editor::parse_path(&out_path).unwrap();
    assert_eq!(names(&f), ["System", "Protein", "SOL", "System_&_!Protein"]);
    assert_eq!(f.groups[3].atoms, (21..=50).collect::<Vec<u32>>());
}

#[test]
fn q_bang_writes_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.ndx");

    let (out, code) = script("small.ndx", "0 & !1\nq!\n");
    assert_eq!(code, 0);
    assert!(out.contains("discarded unsaved changes"));
    assert!(!out_path.exists());
}

#[test]
fn the_input_file_is_never_written_to() {
    let f = fixture("small.ndx");
    let before = bytes(&f);

    // `q` pointed straight back at the input must be refused.
    let (out, code) = script("small.ndx", &format!("0 & !1\nq {}\n", f.display()));
    assert!(out.contains("refusing to overwrite"), "{out}");
    assert_eq!(bytes(&f), before);
    // The session then hits EOF while still dirty and cannot save, so it exits non-zero.
    assert_eq!(code, 3);
}

#[test]
fn q_without_a_destination_fails_loudly_rather_than_silently_dropping_work() {
    let (out, code) = script("small.ndx", "0 & !1\nq\n");
    assert!(out.contains("no output file set"), "{out}");
    assert!(out.contains("q FILE"));
    assert_eq!(code, 3);
}

#[test]
fn dash_o_gives_q_a_destination() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.ndx");

    ndxed()
        .args([
            "edit",
            fixture("small.ndx").to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
        ])
        .write_stdin("0 & !1\nq\n")
        .assert()
        .success()
        .stdout(predicate::str::contains("wrote 4 group(s)"));

    assert_eq!(ndx_editor::parse_path(&out_path).unwrap().len(), 4);
}

/// EOF (Ctrl-D, or the end of a piped script) behaves like `q`.
#[test]
fn eof_saves_like_q() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.ndx");

    ndxed()
        .args([
            "edit",
            fixture("small.ndx").to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
        ])
        .write_stdin("0 & !1\n") // no `q`
        .assert()
        .success();

    assert_eq!(ndx_editor::parse_path(&out_path).unwrap().len(), 4);
}

#[test]
fn eof_with_no_changes_writes_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.ndx");

    ndxed()
        .args([
            "edit",
            fixture("small.ndx").to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
        ])
        .write_stdin("l\n")
        .assert()
        .success();

    assert!(!out_path.exists(), "a read-only session must not write");
}

#[test]
fn undo_and_redo() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.ndx");

    let (out, _) = script(
        "small.ndx",
        &format!("0 & !1\nname 3 Solvent\nundo\nundo\nq {}\n", out_path.display()),
    );
    assert!(out.contains("Solvent"), "{out}");

    // Both edits undone: back to the three original groups.
    let f = ndx_editor::parse_path(&out_path).unwrap();
    assert_eq!(names(&f), ["System", "Protein", "SOL"]);

    let (out, _) = script(
        "small.ndx",
        &format!("0 & !1\nundo\nredo\nq {}\n", out_path.display()),
    );
    assert!(!out.contains("error"), "{out}");
    assert_eq!(ndx_editor::parse_path(&out_path).unwrap().len(), 4);
}

/// A command that fails must leave no snapshot behind, or the user would have to press `undo`
/// once to unwind the failure before reaching their last real edit.
#[test]
fn a_failed_command_does_not_consume_an_undo() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.ndx");

    let (out, code) = script(
        "small.ndx",
        &format!("0 & !1\nNoSuchGroup\nundo\nq {}\n", out_path.display()),
    );
    assert!(out.contains("no group named"), "{out}");
    assert_eq!(code, 0, "a bad command mid-session is not fatal");

    // The single `undo` must have unwound the *selection*, not the failure.
    let f = ndx_editor::parse_path(&out_path).unwrap();
    assert_eq!(f.len(), 3);
}

#[test]
fn rename_del_keep_head_and_split() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.ndx");

    let (out, code) = script(
        "small.ndx",
        &format!(
            "name 1 Prot\nhead 1 5 First5\nsplit 2 parts 2\ndel 0\nq {}\n",
            out_path.display()
        ),
    );
    assert_eq!(code, 0, "{out}");

    let f = ndx_editor::parse_path(&out_path).unwrap();
    assert_eq!(names(&f), ["Prot", "SOL", "First5", "SOL_1", "SOL_2"]);
    assert_eq!(f.groups[2].atoms, [1, 2, 3, 4, 5]);
    assert_eq!(f.groups[3].len() + f.groups[4].len(), 30);
}

#[test]
fn atoms_and_keep() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.ndx");

    let (_, code) = script(
        "small.ndx",
        &format!("a 1-3,7\nkeep 3\nq {}\n", out_path.display()),
    );
    assert_eq!(code, 0);

    let f = ndx_editor::parse_path(&out_path).unwrap();
    assert_eq!(f.len(), 1);
    assert_eq!(f.groups[0].atoms, [1, 2, 3, 7]);
}

#[test]
fn merge_pulls_in_another_file() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.ndx");

    let (out, code) = script(
        "small.ndx",
        &format!(
            "merge {}\nq {}\n",
            fixture("other.ndx").display(),
            out_path.display()
        ),
    );
    assert_eq!(code, 0, "{out}");
    assert!(out.contains("duplicate group name"), "{out}");

    let f = ndx_editor::parse_path(&out_path).unwrap();
    assert_eq!(names(&f), ["System", "Protein", "SOL", "Protein", "Ligand"]);
}

#[test]
fn diff_against_a_file() {
    let (out, code) = script(
        "small.ndx",
        &format!("diff {}\nq!\n", fixture("other.ndx").display()),
    );
    assert_eq!(code, 0);
    assert!(out.contains("- System"), "{out}");
    assert!(out.contains("+ Ligand"));
}

/// Enter on an empty line re-prints the table — the thing you want most often after an edit.
#[test]
fn a_bare_enter_lists_the_groups() {
    let (out, code) = script("small.ndx", "0 & !1\n\nq!\n");
    assert_eq!(code, 0);

    // The table is printed twice: once in the opening banner, once for the blank line.
    assert_eq!(out.matches("  0 System").count(), 2, "{out}");

    // And the second one includes the group the expression had just added.
    let after_blank = out.split("> \n").nth(1).expect("a blank-line prompt");
    assert!(after_blank.contains("System_&_!Protein"), "{after_blank}");
}

/// A comment stays a no-op, so a piped script can annotate itself without printing a table after
/// every remark.
#[test]
fn a_comment_line_prints_nothing() {
    let (out, _) = script("small.ndx", "# just a note\nq!\n");
    assert_eq!(out.matches("System ").count(), 1, "only the banner\n{out}");
}

#[test]
fn help_lists_the_commands() {
    let (out, _) = script("small.ndx", "help\nq!\n");
    assert!(out.contains("<expression>"));
    assert!(out.contains("split"));
    assert!(out.contains("q!"));
}

#[test]
fn expr_escapes_a_command_word_collision() {
    // `del` here is an expression, not the delete command.
    let (out, _) = script("small.ndx", "expr 1 | 2\nq!\n");
    assert!(out.contains("Protein_|_SOL"), "{out}");
}

#[test]
fn a_structure_expression_explains_itself_without_killing_the_session() {
    let (out, code) = script("small.ndx", "element H & bonded Protein\nl\nq!\n");
    assert!(out.contains("needs a structure file"), "{out}");
    assert!(out.contains("^"), "the caret should point at the expression");
    // The session carried on, and nothing was added.
    assert!(out.contains("SOL"));
    assert_eq!(code, 0);
}

#[test]
fn dry_run_writes_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.ndx");

    ndxed()
        .args([
            "edit",
            fixture("small.ndx").to_str().unwrap(),
            "--dry-run",
            "-o",
            out_path.to_str().unwrap(),
        ])
        .write_stdin("0 & !1\nq\n")
        .assert()
        .success()
        .stdout(predicate::str::contains("would write"));

    assert!(!out_path.exists());
}

/// `ndxed FILE` with no subcommand is the editor.
#[test]
fn a_bare_file_argument_opens_the_editor() {
    ndxed()
        .arg(fixture("small.ndx"))
        .write_stdin("l\nq!\n")
        .assert()
        .success()
        .stdout(predicate::str::contains("Reading index file"));
}
