pub mod command;
pub mod exec;
pub mod session;

use std::io::{BufRead, IsTerminal, Write};
use std::path::PathBuf;

use ndx_editor::error::{NdxError, Result, exit, render_caret};
use ndx_editor::parse::{ParseOptions, parse_path_with};
use ndx_editor::universe::UniverseSpec;

use crate::cli::render;
use crate::repl::command::parse_line;
use crate::repl::exec::{Flow, exec};
use crate::repl::session::Session;

pub struct Options {
    pub file: PathBuf,
    pub out: Option<PathBuf>,
    pub force_overwrite: bool,
    pub dry_run: bool,
    pub no_readline: bool,
    pub universe: UniverseSpec,
    pub quiet: bool,
}

pub fn run(opts: Options) -> Result<i32> {
    let (ndx, report) = parse_path_with(&opts.file, &ParseOptions::default())?;

    let mut out = anstream::stdout();
    writeln!(out, "Reading index file '{}'", opts.file.display()).ok();
    if !opts.quiet {
        for name in &report.unsorted {
            writeln!(out, "note: group {name:?} is not sorted (kept as-is)").ok();
        }
        for name in &report.duplicated {
            writeln!(out, "note: group {name:?} has duplicate atoms (kept as-is)").ok();
        }
    }
    writeln!(out).ok();
    render::list(&mut out, &ndx, false, false).ok();
    writeln!(out, "\nType `help` for commands, `q FILE` to save and quit.").ok();
    if opts.out.is_none() {
        writeln!(
            out,
            "No output file set: `q` will ask for one (or use `q FILE`)."
        )
        .ok();
    }
    let _ = out.flush();

    let mut session = Session::open(
        opts.file.clone(),
        ndx,
        opts.out,
        opts.universe,
        opts.dry_run,
        opts.force_overwrite,
    );

    let interactive = std::io::stdin().is_terminal() && !opts.no_readline;
    if interactive {
        #[cfg(feature = "repl")]
        {
            return readline_loop(&mut session);
        }
    }
    piped_loop(&mut session)
}

/// The driver used when stdin is not a terminal (and when built without the `repl` feature).
///
/// Keeping this the primary implementation is what makes the REPL testable without a pty:
/// `printf 'l\n0 & !1\nq out.ndx\n' | ndxed edit in.ndx` exercises the real code path.
fn piped_loop(session: &mut Session) -> Result<i32> {
    let stdin = std::io::stdin();
    let mut out = anstream::stdout();

    for line in stdin.lock().lines() {
        let line = line.map_err(|e| NdxError::io("<stdin>", e))?;
        writeln!(out, "> {line}").ok();
        match step(session, &line, &mut out) {
            Ok(Flow::Quit) => {
                let _ = out.flush();
                return Ok(exit::OK);
            }
            Ok(Flow::Continue) => {}
            Err(e) => {
                // A bad command is not fatal mid-session; report it and keep going.
                print_error(&mut out, &line, &e);
            }
        }
        let _ = out.flush();
    }

    // EOF (Ctrl-D, or the end of a piped script) means `q`: save if we can, otherwise say why not.
    finish_at_eof(session, &mut out)
}

fn finish_at_eof(session: &mut Session, out: &mut impl Write) -> Result<i32> {
    if !session.dirty {
        return Ok(exit::OK);
    }
    match session.save(None) {
        Ok(to) => {
            let verb = if session.dry_run {
                "would write"
            } else {
                "wrote"
            };
            writeln!(out, "{verb} {} group(s) to {}", session.ndx.len(), to.display()).ok();
            let _ = out.flush();
            Ok(exit::OK)
        }
        Err(e) => {
            let _ = out.flush();
            Err(e)
        }
    }
}

fn step(session: &mut Session, line: &str, out: &mut impl Write) -> Result<Flow> {
    match parse_line(line)? {
        Some(cmd) => exec(session, cmd, out),
        None => Ok(Flow::Continue),
    }
}

fn print_error(w: &mut impl Write, line: &str, e: &NdxError) {
    writeln!(w, "error: {e}").ok();
    // An expression error knows which byte it choked on; the command word is not part of the
    // expression when `expr ...` was used, so re-find the source in the line.
    if let Some(span) = e.span() {
        let src = line
            .trim()
            .strip_prefix("expr ")
            .map(str::trim)
            .unwrap_or_else(|| line.trim());
        if span.end <= src.len() {
            writeln!(w, "{}", render_caret(src, &span)).ok();
        }
    }
}

#[cfg(feature = "repl")]
fn readline_loop(session: &mut Session) -> Result<i32> {
    use rustyline::error::ReadlineError;

    let mut rl: rustyline::DefaultEditor = rustyline::DefaultEditor::new()
        .map_err(|e| NdxError::Other(format!("could not start the line editor: {e}")))?;
    let history = dirs_history();
    if let Some(h) = &history {
        let _ = rl.load_history(h);
    }

    let mut out = anstream::stdout();
    let mut interrupts = 0u8;

    loop {
        match rl.readline("> ") {
            Ok(line) => {
                interrupts = 0;
                let _ = rl.add_history_entry(line.as_str());
                match step(session, &line, &mut out) {
                    Ok(Flow::Quit) => break,
                    Ok(Flow::Continue) => {}
                    Err(e) => print_error(&mut out, &line, &e),
                }
                let _ = out.flush();
            }

            // Ctrl-C cancels the line; twice in a row on an empty line leaves without saving.
            Err(ReadlineError::Interrupted) => {
                interrupts += 1;
                if interrupts >= 2 {
                    writeln!(out, "quitting without saving").ok();
                    break;
                }
                writeln!(out, "(interrupted; press Ctrl-C again to quit without saving)").ok();
            }

            // Ctrl-D is `q`.
            Err(ReadlineError::Eof) => {
                if let Err(e) = finish_at_eof(session, &mut out) {
                    writeln!(out, "error: {e}").ok();
                    let _ = out.flush();
                    if let Some(to) = ask_where_to_save(&mut rl, session, &mut out)? {
                        writeln!(out, "wrote {} group(s) to {}", session.ndx.len(), to.display())
                            .ok();
                    }
                }
                break;
            }

            Err(e) => return Err(NdxError::Other(format!("line editor: {e}"))),
        }
    }

    if let Some(h) = &history {
        let _ = rl.save_history(h);
    }
    let _ = out.flush();
    Ok(exit::OK)
}

/// `q` with no destination and no `-o`: ask, rather than silently dropping the work.
#[cfg(feature = "repl")]
fn ask_where_to_save(
    rl: &mut rustyline::DefaultEditor,
    session: &mut Session,
    out: &mut impl Write,
) -> Result<Option<PathBuf>> {
    use rustyline::error::ReadlineError;

    loop {
        match rl.readline("Save to (empty to discard): ") {
            Ok(line) => {
                let line = line.trim();
                if line.is_empty() {
                    writeln!(out, "discarded unsaved changes").ok();
                    return Ok(None);
                }
                match session.save(Some(std::path::Path::new(line))) {
                    Ok(to) => return Ok(Some(to)),
                    Err(e) => writeln!(out, "error: {e}").ok(),
                };
            }
            Err(ReadlineError::Interrupted | ReadlineError::Eof) => {
                writeln!(out, "discarded unsaved changes").ok();
                return Ok(None);
            }
            Err(e) => return Err(NdxError::Other(format!("line editor: {e}"))),
        }
    }
}

#[cfg(feature = "repl")]
fn dirs_history() -> Option<PathBuf> {
    std::env::var_os("HOME").map(|h| PathBuf::from(h).join(".ndxed_history"))
}
