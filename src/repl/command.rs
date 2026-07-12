use std::path::PathBuf;

use ndx_editor::error::{NdxError, Result};
use ndx_editor::model::GroupRef;
use ndx_editor::ops::End;
use ndx_editor::ops::split::How;

#[derive(Clone, Debug, PartialEq)]
pub enum Cmd {
    List { long: bool },
    Name(GroupRef, String),
    Del(Vec<GroupRef>),
    Keep(Vec<GroupRef>),
    Atoms(String),
    Head(GroupRef, usize, End, Option<String>),
    Split(GroupRef, How, Option<String>, bool),
    Merge(PathBuf),
    Diff(PathBuf),
    Undo,
    Redo,
    Write(Option<PathBuf>),
    Quit(Option<PathBuf>),
    QuitNoSave,
    Help(Option<String>),
    /// A bare expression: `0 & !1`, `Protein | SOL`, ...
    Expr(String),
}

/// Every word that is a command rather than the start of an expression.
///
/// The dispatch rule is deliberately the same as make_ndx's: if the first token is one of these,
/// the line is a command; otherwise the whole line is an expression, so `0 & !1` works unchanged.
/// `expr <...>` is the escape hatch for a group whose name collides with a keyword.
const KEYWORDS: &[&str] = &[
    "l", "list", "name", "del", "delete", "keep", "a", "atoms", "head", "tail", "split", "merge",
    "diff", "undo", "redo", "w", "write", "q", "quit", "q!", "h", "help", "?", "expr",
];

pub fn parse_line(line: &str) -> Result<Option<Cmd>> {
    let line = line.trim();
    if line.is_empty() || line.starts_with('#') {
        return Ok(None);
    }

    let (head, rest) = match line.split_once(char::is_whitespace) {
        Some((h, r)) => (h, r.trim()),
        None => (line, ""),
    };

    if !KEYWORDS.contains(&head) {
        return Ok(Some(Cmd::Expr(line.to_string())));
    }

    let cmd = match head {
        "l" | "list" => Cmd::List {
            long: rest == "--long" || rest == "-l",
        },

        "name" => {
            let (g, new) = rest
                .split_once(char::is_whitespace)
                .ok_or_else(|| usage("name <group> <new name>"))?;
            let new = new.trim();
            if new.is_empty() {
                return Err(usage("name <group> <new name>"));
            }
            Cmd::Name(g.parse()?, new.to_string())
        }

        "del" | "delete" => Cmd::Del(refs(rest, "del <group> [<group>...]")?),
        "keep" => Cmd::Keep(refs(rest, "keep <group> [<group>...]")?),

        "a" | "atoms" => {
            if rest.is_empty() {
                return Err(usage("a <ranges>, e.g. `a 1-10,15,20-30`"));
            }
            Cmd::Atoms(rest.to_string())
        }

        "head" | "tail" => {
            let end = if head == "head" { End::Head } else { End::Tail };
            let mut it = rest.split_ascii_whitespace();
            let g = it
                .next()
                .ok_or_else(|| usage("head <group> <N> [name]"))?
                .parse()?;
            let n: usize = it
                .next()
                .ok_or_else(|| usage("head <group> <N> [name]"))?
                .parse()
                .map_err(|_| usage("head <group> <N> [name]"))?;
            let name: Vec<&str> = it.collect();
            let name = (!name.is_empty()).then(|| name.join(" "));
            Cmd::Head(g, n, end, name)
        }

        "split" => parse_split(rest)?,

        "merge" => {
            if rest.is_empty() {
                return Err(usage("merge <file.ndx>"));
            }
            Cmd::Merge(PathBuf::from(rest))
        }

        "diff" => {
            if rest.is_empty() {
                return Err(usage("diff <file.ndx>"));
            }
            Cmd::Diff(PathBuf::from(rest))
        }

        "undo" => Cmd::Undo,
        "redo" => Cmd::Redo,

        "w" | "write" => Cmd::Write(path_or_none(rest)),
        "q" | "quit" => Cmd::Quit(path_or_none(rest)),
        "q!" => Cmd::QuitNoSave,

        "h" | "help" | "?" => Cmd::Help((!rest.is_empty()).then(|| rest.to_string())),

        "expr" => {
            if rest.is_empty() {
                return Err(usage("expr <expression>"));
            }
            Cmd::Expr(rest.to_string())
        }

        _ => unreachable!("head was checked against KEYWORDS"),
    };

    Ok(Some(cmd))
}

fn parse_split(rest: &str) -> Result<Cmd> {
    const USAGE: &str = "split <group> parts|size|at <N> [prefix] [--replace]";
    let mut words: Vec<&str> = rest.split_ascii_whitespace().collect();

    let replace = words.contains(&"--replace");
    words.retain(|w| *w != "--replace");

    if words.len() < 3 {
        return Err(usage(USAGE));
    }
    let g: GroupRef = words[0].parse()?;
    let how = match words[1] {
        "parts" => How::Parts(words[2].parse().map_err(|_| usage(USAGE))?),
        "size" => How::Size(words[2].parse().map_err(|_| usage(USAGE))?),
        "at" => {
            let mut cuts = Vec::new();
            for tok in words[2].split(',').filter(|t| !t.is_empty()) {
                cuts.push(tok.parse().map_err(|_| usage(USAGE))?);
            }
            if cuts.is_empty() {
                return Err(usage(USAGE));
            }
            How::At(cuts)
        }
        other => {
            return Err(NdxError::Other(format!(
                "unknown split mode {other:?}\nusage: {USAGE}"
            )));
        }
    };
    let prefix = words.get(3).map(|s| (*s).to_string());
    Ok(Cmd::Split(g, how, prefix, replace))
}

fn refs(rest: &str, usage_str: &str) -> Result<Vec<GroupRef>> {
    let out: Vec<GroupRef> = rest
        .split_ascii_whitespace()
        .map(str::parse)
        .collect::<Result<_>>()?;
    if out.is_empty() {
        return Err(usage(usage_str));
    }
    Ok(out)
}

fn path_or_none(rest: &str) -> Option<PathBuf> {
    (!rest.is_empty()).then(|| PathBuf::from(rest))
}

fn usage(s: &str) -> NdxError {
    NdxError::Other(format!("usage: {s}"))
}

pub const HELP: &str = "\
  l | list [--long]           list the groups
  <expression>                add a group, e.g. `0 & !1`, `Protein | SOL`, `(1|2) \\ 3`
  expr <expression>           same, for a group whose name collides with a command
  a | atoms <ranges>          add a group from atom numbers, e.g. `a 1-10,15,20-30`
  name <group> <new name>     rename a group
  del <group>...              delete groups (ids are renumbered)
  keep <group>...             keep only these groups, in the order given
  head <group> <N> [name]     add the first N atoms of a group (in file order)
  tail <group> <N> [name]     add the last N atoms
  split <group> parts <K>     split into K near-equal parts
  split <group> size <K>      split into chunks of K atoms
  split <group> at 100,250    cut after those positions
                              ... all three take [prefix] [--replace]
  merge <file.ndx>            append another file's groups to this session
  diff <file.ndx>             compare this session against a file
  undo | redo                 step through the edit history
  w [file]                    write now (and remember the destination)
  q [file]                    write and quit
  q!                          quit without writing
  h | help                    this
\n\
Operators, tightest first:  !  &  \\  |
Groups: by id (`3`), name (`SOL`), duplicate (`SOL#2`), or wildcard (`Fiber*`, `O?`).
A wildcard means every group whose name matches — in an expression, their union.
Quote it to match a name literally: `\"O*\"`.
The input file is never overwritten: `q` needs -o or a filename.";

#[cfg(test)]
mod tests {
    use super::*;

    fn p(line: &str) -> Cmd {
        parse_line(line).unwrap().unwrap()
    }

    #[test]
    fn blank_and_comment_lines_do_nothing() {
        assert!(parse_line("").unwrap().is_none());
        assert!(parse_line("   ").unwrap().is_none());
        assert!(parse_line("# note").unwrap().is_none());
    }

    /// The whole point of the dispatch rule: make_ndx muscle memory keeps working.
    #[test]
    fn a_bare_expression_stays_an_expression() {
        assert_eq!(p("0 & !1"), Cmd::Expr("0 & !1".into()));
        assert_eq!(p("Protein | SOL"), Cmd::Expr("Protein | SOL".into()));
    }

    #[test]
    fn list() {
        assert_eq!(p("l"), Cmd::List { long: false });
        assert_eq!(p("list --long"), Cmd::List { long: true });
    }

    #[test]
    fn rename_takes_a_name_with_spaces() {
        assert_eq!(
            p("name 3 C alpha"),
            Cmd::Name(GroupRef::Id(3), "C alpha".into())
        );
    }

    #[test]
    fn del_and_keep_take_several() {
        assert_eq!(p("del 1 2"), Cmd::Del(vec![GroupRef::Id(1), GroupRef::Id(2)]));
        assert_eq!(
            p("keep 0 SOL"),
            Cmd::Keep(vec![GroupRef::Id(0), GroupRef::name("SOL")])
        );
    }

    #[test]
    fn head_and_tail() {
        assert_eq!(
            p("head 1 100"),
            Cmd::Head(GroupRef::Id(1), 100, End::Head, None)
        );
        assert_eq!(
            p("tail 1 5 Last"),
            Cmd::Head(GroupRef::Id(1), 5, End::Tail, Some("Last".into()))
        );
    }

    #[test]
    fn split_modes() {
        assert_eq!(
            p("split 1 parts 4"),
            Cmd::Split(GroupRef::Id(1), How::Parts(4), None, false)
        );
        assert_eq!(
            p("split 1 size 10 chunk --replace"),
            Cmd::Split(GroupRef::Id(1), How::Size(10), Some("chunk".into()), true)
        );
        assert_eq!(
            p("split 1 at 100,250"),
            Cmd::Split(GroupRef::Id(1), How::At(vec![100, 250]), None, false)
        );
    }

    #[test]
    fn quit_and_write_take_an_optional_path() {
        assert_eq!(p("q"), Cmd::Quit(None));
        assert_eq!(p("q out.ndx"), Cmd::Quit(Some("out.ndx".into())));
        assert_eq!(p("q!"), Cmd::QuitNoSave);
        assert_eq!(p("w"), Cmd::Write(None));
        assert_eq!(p("w out.ndx"), Cmd::Write(Some("out.ndx".into())));
    }

    #[test]
    fn expr_escapes_a_keyword_collision() {
        // A group actually called `del` can still be selected.
        assert_eq!(p("expr del | 1"), Cmd::Expr("del | 1".into()));
    }

    #[test]
    fn bad_usage_is_reported() {
        assert!(parse_line("name").is_err());
        assert!(parse_line("name 3").is_err());
        assert!(parse_line("head 1").is_err());
        assert!(parse_line("head 1 x").is_err());
        assert!(parse_line("split 1 parts").is_err());
        assert!(parse_line("split 1 bogus 4").is_err());
        assert!(parse_line("del").is_err());
        assert!(parse_line("merge").is_err());
    }
}
