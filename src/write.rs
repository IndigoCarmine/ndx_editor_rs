use std::io::Write;
use std::path::{Path, PathBuf};

use crate::error::{NdxError, Result};
use crate::model::IndexFile;

#[derive(Clone, Copy, Debug)]
pub struct WriteOptions {
    pub per_line: usize,
    pub width: usize,
}

impl Default for WriteOptions {
    fn default() -> Self {
        // Matches `write_index()` in GROMACS: 15 indices per line, "%4d".
        WriteOptions {
            per_line: 15,
            width: 4,
        }
    }
}

/// Write an index file in GROMACS' exact byte format.
///
/// GROMACS emits `[ name ]`, then for each atom a separator (`\n` every 15th, else a space)
/// followed by `%4d`, then a final `\n`. Reproducing that exactly is what makes
/// `write(parse(f)) == f` hold for gmx-produced files.
pub fn write<W: Write>(w: &mut W, ndx: &IndexFile, o: &WriteOptions) -> std::io::Result<()> {
    let per_line = o.per_line.max(1);
    for g in &ndx.groups {
        write!(w, "[ {} ]", g.name)?;
        for (k, atom) in g.atoms.iter().enumerate() {
            let sep = if k % per_line == 0 { '\n' } else { ' ' };
            write!(w, "{sep}{atom:>width$}", width = o.width)?;
        }
        writeln!(w)?;
    }
    Ok(())
}

pub fn to_string(ndx: &IndexFile, o: &WriteOptions) -> String {
    let mut buf: Vec<u8> = Vec::new();
    write(&mut buf, ndx, o).expect("writing to a Vec cannot fail");
    String::from_utf8(buf).expect("group names and indices are valid UTF-8")
}

/// Where a mutating command sends its result.
///
/// The input file is never a valid target: every command defaults to stdout, and `-o` names a
/// different file. `Sink::Path` refuses to clobber `guard` unless `force` was passed.
#[derive(Clone, Debug)]
pub enum Sink {
    Stdout,
    Path(PathBuf),
}

impl Sink {
    pub fn new(out: Option<PathBuf>) -> Self {
        match out {
            Some(p) if p != Path::new("-") => Sink::Path(p),
            _ => Sink::Stdout,
        }
    }

    /// Write `ndx` out. `guards` are the input paths we must not overwrite.
    pub fn commit(
        &self,
        ndx: &IndexFile,
        o: &WriteOptions,
        guards: &[&Path],
        force: bool,
    ) -> Result<()> {
        match self {
            Sink::Stdout => {
                let stdout = std::io::stdout();
                let mut lock = stdout.lock();
                write(&mut lock, ndx, o).map_err(|e| NdxError::io("<stdout>", e))?;
                lock.flush().map_err(|e| NdxError::io("<stdout>", e))
            }
            Sink::Path(p) => {
                if !force && guards.iter().any(|g| same_file(p, g)) {
                    return Err(NdxError::WouldOverwriteInput { path: p.clone() });
                }
                write_atomic(p, ndx, o)
            }
        }
    }
}

/// Two paths point at the same file. Falls back to a textual compare when either side can't be
/// canonicalized (e.g. the output does not exist yet — in which case it is not the input).
fn same_file(a: &Path, b: &Path) -> bool {
    match (a.canonicalize(), b.canonicalize()) {
        (Ok(a), Ok(b)) => a == b,
        _ => a == b,
    }
}

/// Write via a temp file in the same directory, then rename. A crash or a mid-write error
/// therefore never leaves a truncated index file behind.
pub fn write_atomic(path: &Path, ndx: &IndexFile, o: &WriteOptions) -> Result<()> {
    let dir = path.parent().filter(|p| !p.as_os_str().is_empty());
    let dir = dir.unwrap_or(Path::new("."));

    let mut tmp = tempfile::NamedTempFile::new_in(dir).map_err(|e| NdxError::io(dir, e))?;
    write(&mut tmp, ndx, o).map_err(|e| NdxError::io(path, e))?;
    tmp.flush().map_err(|e| NdxError::io(path, e))?;
    tmp.persist(path)
        .map_err(|e| NdxError::io(path, e.error))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Group;
    use crate::parse::parse_str;

    #[test]
    fn gmx_byte_format() {
        let f = IndexFile {
            groups: vec![Group::new("System", (1..=17).collect())],
        };
        let s = to_string(&f, &WriteOptions::default());
        let expected = "[ System ]\n\
             \x20  1    2    3    4    5    6    7    8    9   10   11   12   13   14   15\n\
             \x20 16   17\n";
        assert_eq!(s, expected);
    }

    #[test]
    fn empty_group_is_just_a_header() {
        let f = IndexFile {
            groups: vec![Group::new("Empty", vec![])],
        };
        assert_eq!(to_string(&f, &WriteOptions::default()), "[ Empty ]\n");
    }

    #[test]
    fn wide_indices_are_not_truncated() {
        let f = IndexFile {
            groups: vec![Group::new("Big", vec![1, 123456])],
        };
        assert_eq!(to_string(&f, &WriteOptions::default()), "[ Big ]\n   1 123456\n");
    }

    #[test]
    fn order_and_duplicates_survive() {
        let f = IndexFile {
            groups: vec![Group::new("Odd", vec![3, 1, 1])],
        };
        let s = to_string(&f, &WriteOptions::default());
        assert_eq!(parse_str(&s, "t").unwrap(), f);
    }

    #[test]
    fn roundtrip_through_parse() {
        let src = "[ System ]\n\
             \x20  1    2    3    4    5    6    7    8    9   10   11   12   13   14   15\n\
             \x20 16   17\n\
             [ C-alpha ]\n\
             \x20  2    5\n";
        let f = parse_str(src, "t").unwrap();
        // The strong one: a gmx-formatted file survives a read/write cycle byte for byte.
        assert_eq!(to_string(&f, &WriteOptions::default()), src);
    }
}
