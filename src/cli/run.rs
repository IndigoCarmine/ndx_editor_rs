use std::io::Write;
use std::path::{Path, PathBuf};

use clap::CommandFactory;

use ndx_editor::error::{NdxError, Result, exit, render_caret};
use ndx_editor::model::IndexFile;
use ndx_editor::ops::diff::diff as make_diff;
use ndx_editor::ops::merge::merge;
use ndx_editor::ops::split::How;
use ndx_editor::ops::{self, End};
use ndx_editor::parse::{ParseOptions, ParseReport, parse_path_with};
use ndx_editor::system::SystemCtx;
use ndx_editor::universe::UniverseSpec;
use ndx_editor::write::{Sink, WriteOptions};

use crate::cli::{Cli, Cmd, OutOpts, UniverseOpts};
use crate::cli::render;
use crate::repl;

/// stderr is for diagnostics; stdout is data, always.
pub struct Ui {
    pub quiet: bool,
}

impl Ui {
    pub fn warn(&self, msg: &str) {
        if !self.quiet {
            let mut err = anstream::stderr();
            let _ = writeln!(err, "warning: {msg}");
        }
    }

    pub fn note(&self, msg: &str) {
        if !self.quiet {
            let mut err = anstream::stderr();
            let _ = writeln!(err, "{msg}");
        }
    }

    fn report(&self, r: &ParseReport, origin: &str) {
        if !r.unsorted.is_empty() {
            self.warn(&format!(
                "{origin}: group(s) with unsorted atoms: {} \
                 (kept as-is; `ndxed fmt --sort` canonicalizes)",
                r.unsorted.join(", ")
            ));
        }
        if !r.duplicated.is_empty() {
            self.warn(&format!(
                "{origin}: group(s) with duplicate atoms: {} \
                 (kept as-is; `ndxed fmt --dedup` canonicalizes)",
                r.duplicated.join(", ")
            ));
        }
    }
}

impl OutOpts {
    fn write_options(&self) -> WriteOptions {
        WriteOptions {
            per_line: self.per_line,
            width: self.width,
        }
    }

    fn commit(&self, ndx: &IndexFile, inputs: &[&Path]) -> Result<()> {
        let sink = Sink::new(self.out.clone());
        sink.commit(ndx, &self.write_options(), inputs, self.force_overwrite)?;
        if let Sink::Path(p) = &sink {
            let mut err = anstream::stderr();
            let _ = writeln!(err, "wrote {} group(s) to {}", ndx.len(), p.display());
        }
        Ok(())
    }
}

impl UniverseOpts {
    pub fn spec(&self) -> UniverseSpec {
        match (&self.natoms, &self.universe) {
            (Some(n), _) => UniverseSpec::Natoms(*n),
            (_, Some(g)) => UniverseSpec::Group(g.clone()),
            _ => UniverseSpec::Auto,
        }
    }
}

fn read(ui: &Ui, path: &Path) -> Result<IndexFile> {
    let (ndx, report) = parse_path_with(path, &ParseOptions::default())?;
    ui.report(&report, &path.display().to_string());
    Ok(ndx)
}

/// Returns the process exit code.
pub fn run(cli: Cli) -> i32 {
    let ui = Ui { quiet: cli.quiet };

    // `ndxed FILE` with no subcommand drops straight into the editor.
    let cmd = match (cli.cmd, cli.file) {
        (Some(cmd), _) => cmd,
        (None, Some(file)) => Cmd::Edit {
            file,
            out: None,
            force_overwrite: false,
            dry_run: false,
            no_readline: false,
            universe: UniverseOpts {
                natoms: None,
                universe: None,
            },
        },
        (None, None) => {
            let _ = Cli::command().print_help();
            return exit::OK;
        }
    };

    match dispatch(&ui, cmd) {
        Ok(code) => code,
        // `ndxed select … | head` closes our stdout early. That is normal shell usage, not a
        // failure: die quietly, the way every other Unix tool does.
        Err(e) if is_broken_pipe(&e) => exit::OK,
        Err(e) => {
            report_error(&e);
            e.exit_code()
        }
    }
}

fn is_broken_pipe(e: &NdxError) -> bool {
    matches!(e, NdxError::Io { source, .. } if source.kind() == std::io::ErrorKind::BrokenPipe)
}

fn report_error(e: &NdxError) {
    let mut err = anstream::stderr();
    let _ = writeln!(err, "error: {e}");
    // Expression errors know where they happened; show the caret if we still have the source.
    if let Some(src) = EXPR_SOURCE.with(|s| s.borrow().clone())
        && let Some(span) = e.span()
    {
        let _ = writeln!(err, "{}", render_caret(&src, &span));
    }
}

thread_local! {
    /// The expression currently being evaluated, so an error can point a caret at it.
    static EXPR_SOURCE: std::cell::RefCell<Option<String>> = const { std::cell::RefCell::new(None) };
}

pub fn set_expr_source(src: Option<String>) {
    EXPR_SOURCE.with(|s| *s.borrow_mut() = src);
}

fn dispatch(ui: &Ui, cmd: Cmd) -> Result<i32> {
    match cmd {
        Cmd::List {
            file,
            long,
            ranges,
            format,
        } => {
            let ndx = read(ui, &file)?;
            let mut out = anstream::stdout();
            render::list_any(&mut out, &ndx, format, long, ranges)
                .map_err(|e| NdxError::io("<stdout>", e))?;
            let _ = out.flush();
            Ok(exit::OK)
        }

        Cmd::Rename {
            file,
            group,
            new_name,
            out,
        } => {
            let mut ndx = read(ui, &file)?;
            let old = ops::rename(&mut ndx, &group, &new_name)?;
            ui.note(&format!("renamed {old:?} -> {new_name:?}"));
            out.commit(&ndx, &[&file])?;
            Ok(exit::OK)
        }

        Cmd::Select {
            file,
            expr,
            name,
            keep_only,
            error_on_empty,
            universe,
            out,
        } => {
            let mut ndx = read(ui, &file)?;
            set_expr_source(Some(expr.clone()));
            let sel = ops::select::select(
                &mut ndx,
                &expr,
                name.as_deref(),
                &universe.spec(),
                &SystemCtx::default(),
                error_on_empty,
            )?;
            set_expr_source(None);
            for w in &sel.warnings {
                ui.warn(w);
            }
            ui.note(&format!("{:>3} {} : {} atoms", sel.id, sel.name, sel.len));
            if keep_only {
                ndx.retain_ids(&[sel.id]);
            }
            out.commit(&ndx, &[&file])?;
            Ok(exit::OK)
        }

        Cmd::Atoms {
            file,
            ranges,
            name,
            out,
        } => {
            let mut ndx = read(ui, &file)?;
            let id = ops::atoms(&mut ndx, &ranges, name.as_deref())?;
            let g = ndx.get(id)?;
            ui.note(&format!("{id:>3} {} : {} atoms", g.name, g.len()));
            out.commit(&ndx, &[&file])?;
            Ok(exit::OK)
        }

        Cmd::Del { file, groups, out } => {
            let mut ndx = read(ui, &file)?;
            let names = ops::delete(&mut ndx, &groups)?;
            ui.note(&format!("deleted: {}", names.join(", ")));
            out.commit(&ndx, &[&file])?;
            Ok(exit::OK)
        }

        Cmd::Keep { file, groups, out } => {
            let mut ndx = read(ui, &file)?;
            ops::keep(&mut ndx, &groups)?;
            out.commit(&ndx, &[&file])?;
            Ok(exit::OK)
        }

        Cmd::Head {
            file,
            group,
            n,
            tail,
            name,
            strict,
            out,
        } => {
            let mut ndx = read(ui, &file)?;
            let end = if tail { End::Tail } else { End::Head };
            let id = ops::head(&mut ndx, &group, n, end, name.as_deref(), strict)?;
            let g = ndx.get(id)?;
            ui.note(&format!("{id:>3} {} : {} atoms", g.name, g.len()));
            out.commit(&ndx, &[&file])?;
            Ok(exit::OK)
        }

        Cmd::Split {
            file,
            group,
            parts,
            size,
            at,
            prefix,
            replace,
            out,
        } => {
            let how = match (parts, size, at) {
                (Some(k), None, None) => How::Parts(k),
                (None, Some(k), None) => How::Size(k),
                (None, None, Some(cuts)) => How::At(cuts),
                _ => {
                    return Err(NdxError::Other(
                        "split needs exactly one of --parts, --size or --at".into(),
                    ));
                }
            };
            let mut ndx = read(ui, &file)?;
            let r = ops::split::split(&mut ndx, &group, &how, prefix.as_deref(), replace)?;
            for (name, len) in r.names.iter().zip(&r.sizes) {
                ui.note(&format!("{name} : {len} atoms"));
            }
            out.commit(&ndx, &[&file])?;
            Ok(exit::OK)
        }

        Cmd::Merge {
            files,
            on_conflict,
            prefix_file,
            out,
        } => {
            let mut loaded: Vec<(PathBuf, IndexFile)> = Vec::with_capacity(files.len());
            for f in &files {
                loaded.push((f.clone(), read(ui, f)?));
            }
            let (merged, report) = merge(
                loaded.iter().map(|(p, f)| (p.as_path(), f)),
                on_conflict,
                prefix_file,
            );
            if !report.conflicts.is_empty() {
                ui.warn(&format!(
                    "duplicate group name(s) across files: {} (--on-conflict {:?})",
                    report.conflicts.join(", "),
                    on_conflict
                ));
            }
            // Every input is equally "the input": -o must not clobber any of them.
            let guards: Vec<&Path> = files.iter().map(PathBuf::as_path).collect();
            out.commit(&merged, &guards)?;
            Ok(exit::OK)
        }

        Cmd::Diff {
            a,
            b,
            atoms,
            by,
            no_exit_code,
        } => {
            let (fa, fb) = (read(ui, &a)?, read(ui, &b)?);
            let d = make_diff(&fa, &fb, by);
            let mut out = anstream::stdout();
            render::diff(
                &mut out,
                &d,
                &a.display().to_string(),
                &b.display().to_string(),
                atoms,
            )
            .map_err(|e| NdxError::io("<stdout>", e))?;
            let _ = out.flush();
            Ok(if d.is_empty() || no_exit_code {
                exit::OK
            } else {
                exit::DIFFERENT
            })
        }

        Cmd::Fmt {
            file,
            sort,
            dedup,
            out,
        } => {
            let mut ndx = read(ui, &file)?;
            ops::fmt(&mut ndx, sort, dedup);
            out.commit(&ndx, &[&file])?;
            Ok(exit::OK)
        }

        Cmd::Edit {
            file,
            out,
            force_overwrite,
            dry_run,
            no_readline,
            universe,
        } => repl::run(repl::Options {
            file,
            out,
            force_overwrite,
            dry_run,
            no_readline,
            universe: universe.spec(),
            quiet: ui.quiet,
        }),

        Cmd::Completions { shell } => {
            // Generate into a buffer rather than straight to stdout: clap_complete *panics* on a
            // write error, so `ndxed completions bash | head` would blow up on the broken pipe.
            let mut buf: Vec<u8> = Vec::new();
            let mut cmd = Cli::command();
            clap_complete::generate(shell, &mut cmd, "ndxed", &mut buf);

            let mut stdout = std::io::stdout();
            stdout
                .write_all(&buf)
                .and_then(|()| stdout.flush())
                .map_err(|e| NdxError::io("<stdout>", e))?;
            Ok(exit::OK)
        }
    }
}
