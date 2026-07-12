use std::io::Write;
use std::path::PathBuf;

use ndx_editor::error::Result;
use ndx_editor::ops::diff::{Match, diff as make_diff};
use ndx_editor::ops::merge::{OnConflict, merge};
use ndx_editor::ops::{self, select::select};
use ndx_editor::parse::parse_path;

use crate::cli::render;
use crate::cli::run::set_expr_source;
use crate::repl::command::{Cmd, HELP};
use crate::repl::session::Session;

pub enum Flow {
    Continue,
    Quit,
}

/// Run one REPL command. Anything printed goes to `w` (stdout).
pub fn exec(s: &mut Session, cmd: Cmd, w: &mut impl Write) -> Result<Flow> {
    match cmd {
        Cmd::Help(_) => {
            writeln!(w, "{HELP}").ok();
        }

        Cmd::List { long } => {
            render::list(w, &s.ndx, long, false).ok();
        }

        Cmd::Expr(src) => {
            s.snapshot();
            set_expr_source(Some(src.clone()));
            let result = select(
                &mut s.ndx,
                &src,
                None,
                &s.universe,
                &s.system,
                false,
            );
            set_expr_source(None);
            match result {
                Ok(sel) => {
                    for warn in &sel.warnings {
                        writeln!(w, "warning: {warn}").ok();
                    }
                    writeln!(w, "{:>3} {} : {} atoms", sel.id, sel.name, sel.len).ok();
                    s.mark_dirty();
                }
                Err(e) => {
                    // A failed command must not leave a snapshot behind, or `undo` would be a
                    // no-op the user has to press twice.
                    s.rollback();
                    return Err(e);
                }
            }
        }

        Cmd::Atoms(ranges) => {
            s.snapshot();
            match ops::atoms(&mut s.ndx, &ranges, None) {
                Ok(id) => {
                    let g = s.ndx.get(id)?;
                    writeln!(w, "{id:>3} {} : {} atoms", g.name, g.len()).ok();
                    s.mark_dirty();
                }
                Err(e) => {
                    s.rollback();
                    return Err(e);
                }
            }
        }

        Cmd::Name(g, new) => {
            s.snapshot();
            match ops::rename(&mut s.ndx, &g, &new) {
                Ok(old) => {
                    writeln!(w, "renamed {old:?} -> {new:?}").ok();
                    s.mark_dirty();
                }
                Err(e) => {
                    s.rollback();
                    return Err(e);
                }
            }
        }

        Cmd::Del(groups) => {
            s.snapshot();
            match ops::delete(&mut s.ndx, &groups) {
                Ok(names) => {
                    writeln!(w, "deleted: {}", names.join(", ")).ok();
                    s.mark_dirty();
                    render::list(w, &s.ndx, false, false).ok();
                }
                Err(e) => {
                    s.rollback();
                    return Err(e);
                }
            }
        }

        Cmd::Keep(groups) => {
            s.snapshot();
            match ops::keep(&mut s.ndx, &groups) {
                Ok(()) => {
                    s.mark_dirty();
                    render::list(w, &s.ndx, false, false).ok();
                }
                Err(e) => {
                    s.rollback();
                    return Err(e);
                }
            }
        }

        Cmd::Head(g, n, end, name) => {
            s.snapshot();
            match ops::head(&mut s.ndx, &g, n, end, name.as_deref(), false) {
                Ok(id) => {
                    let g = s.ndx.get(id)?;
                    writeln!(w, "{id:>3} {} : {} atoms", g.name, g.len()).ok();
                    s.mark_dirty();
                }
                Err(e) => {
                    s.rollback();
                    return Err(e);
                }
            }
        }

        Cmd::Split(g, how, prefix, replace) => {
            s.snapshot();
            match ops::split::split(&mut s.ndx, &g, &how, prefix.as_deref(), replace, &s.system) {
                Ok(r) => {
                    for (name, len) in r.names.iter().zip(&r.sizes) {
                        writeln!(w, "    {name} : {len} atoms").ok();
                    }
                    s.mark_dirty();
                    render::list(w, &s.ndx, false, false).ok();
                }
                Err(e) => {
                    s.rollback();
                    return Err(e);
                }
            }
        }

        Cmd::Merge(path) => {
            let other = parse_path(&path)?;
            s.snapshot();
            let here = std::mem::take(&mut s.ndx);
            let files = [
                (s.source.as_path(), &here),
                (path.as_path(), &other),
            ];
            let (merged, report) = merge(files, OnConflict::KeepBoth, false);
            s.ndx = merged;
            if !report.conflicts.is_empty() {
                writeln!(
                    w,
                    "warning: duplicate group name(s): {} (use `name <id> ...` to disambiguate)",
                    report.conflicts.join(", ")
                )
                .ok();
            }
            s.mark_dirty();
            render::list(w, &s.ndx, false, false).ok();
        }

        Cmd::Diff(path) => {
            let other = parse_path(&path)?;
            let d = make_diff(&s.ndx, &other, Match::Name);
            render::diff(w, &d, "(session)", &path.display().to_string(), true).ok();
        }

        Cmd::Undo => {
            s.undo()?;
            render::list(w, &s.ndx, false, false).ok();
        }

        Cmd::Redo => {
            s.redo()?;
            render::list(w, &s.ndx, false, false).ok();
        }

        Cmd::Write(path) => {
            let to = save(s, path, w)?;
            let _ = to;
        }

        Cmd::Quit(path) => {
            save(s, path, w)?;
            return Ok(Flow::Quit);
        }

        Cmd::QuitNoSave => {
            if s.dirty {
                writeln!(w, "discarded unsaved changes").ok();
            }
            return Ok(Flow::Quit);
        }
    }

    Ok(Flow::Continue)
}

fn save(s: &mut Session, path: Option<PathBuf>, w: &mut impl Write) -> Result<PathBuf> {
    let to = s.save(path.as_deref())?;
    let verb = if s.dry_run { "would write" } else { "wrote" };
    writeln!(w, "{verb} {} group(s) to {}", s.ndx.len(), to.display()).ok();
    Ok(to)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndx_editor::model::{Group, IndexFile};
    use ndx_editor::universe::UniverseSpec;

    use crate::repl::command::parse_line;

    fn session() -> Session {
        let ndx = IndexFile {
            groups: vec![
                Group::new("System", (1..=10).collect()),
                Group::new("Protein", vec![1, 2, 3]),
                Group::new("SOL", vec![4, 5]),
            ],
        };
        Session::open(
            PathBuf::from("in.ndx"),
            ndx,
            None,
            UniverseSpec::Auto,
            false,
            false,
        )
    }

    fn run(s: &mut Session, line: &str) -> Result<String> {
        let mut buf: Vec<u8> = Vec::new();
        let cmd = parse_line(line)?.expect("a command");
        exec(s, cmd, &mut buf)?;
        Ok(String::from_utf8(buf).unwrap())
    }

    #[test]
    fn an_expression_appends_a_group() {
        let mut s = session();
        let out = run(&mut s, "1 & !2").unwrap();
        assert!(out.contains("Protein_&_!SOL"));
        assert_eq!(s.ndx.len(), 4);
        assert!(s.dirty);
    }

    #[test]
    fn a_failed_command_leaves_no_snapshot_to_undo() {
        let mut s = session();
        assert!(run(&mut s, "NoSuchGroup").is_err());
        assert_eq!(s.ndx.len(), 3);
        // If the failure had left its snapshot behind, this would silently succeed.
        assert!(s.undo().is_err(), "a failed command must not push undo state");
    }

    #[test]
    fn undo_after_a_real_edit_works() {
        let mut s = session();
        run(&mut s, "1 & !2").unwrap();
        run(&mut s, "undo").unwrap();
        assert_eq!(s.ndx.len(), 3);
    }

    #[test]
    fn rename_then_undo() {
        let mut s = session();
        run(&mut s, "name 1 Prot").unwrap();
        assert_eq!(s.ndx.groups[1].name, "Prot");
        run(&mut s, "undo").unwrap();
        assert_eq!(s.ndx.groups[1].name, "Protein");
    }

    #[test]
    fn del_renumbers_and_undo_restores() {
        let mut s = session();
        run(&mut s, "del 0").unwrap();
        assert_eq!(s.ndx.groups[0].name, "Protein");
        run(&mut s, "undo").unwrap();
        assert_eq!(s.ndx.groups[0].name, "System");
    }

    #[test]
    fn quit_saves_and_stops() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.ndx");
        let mut s = session();
        run(&mut s, "1 & !2").unwrap();

        let cmd = parse_line(&format!("q {}", out.display())).unwrap().unwrap();
        let mut buf = Vec::new();
        assert!(matches!(exec(&mut s, cmd, &mut buf).unwrap(), Flow::Quit));

        let written = parse_path(&out).unwrap();
        assert_eq!(written.len(), 4);
        assert_eq!(written.groups[3].name, "Protein_&_!SOL");
    }

    #[test]
    fn quit_without_a_destination_errors() {
        let mut s = session();
        assert!(run(&mut s, "q").is_err());
    }

    #[test]
    fn structure_features_report_what_is_missing() {
        let mut s = session();
        let err = run(&mut s, "element H & bonded 1").unwrap_err();
        assert!(err.to_string().contains(".gro"), "{err}");
        assert_eq!(s.ndx.len(), 3, "nothing should have been added");
    }

    #[test]
    fn head_uses_file_order() {
        let mut s = session();
        s.ndx.push(Group::new("Odd", vec![9, 7, 8]));
        run(&mut s, "head 3 2").unwrap();
        assert_eq!(s.ndx.groups[4].atoms, [9, 7]);
    }

    #[test]
    fn split_lists_the_parts() {
        let mut s = session();
        let out = run(&mut s, "split 0 parts 2").unwrap();
        assert!(out.contains("System_1"));
        assert!(out.contains("System_2"));
    }
}
