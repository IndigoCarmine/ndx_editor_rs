use std::io::Read;
use std::path::Path;

use crate::error::{NdxError, ParseError, Result};
use crate::model::{AtomId, Group, IndexFile};

#[derive(Clone, Debug, Default)]
pub struct ParseOptions {
    /// Treat unsorted groups and duplicate atoms as errors instead of just noting them.
    pub strict: bool,
}

/// What we noticed while reading, so the caller can warn once on stderr.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ParseReport {
    pub unsorted: Vec<String>,
    pub duplicated: Vec<String>,
}

impl ParseReport {
    pub fn is_clean(&self) -> bool {
        self.unsorted.is_empty() && self.duplicated.is_empty()
    }
}

pub fn parse_str(src: &str, origin: &str) -> Result<IndexFile> {
    parse_with(src, origin, &ParseOptions::default()).map(|(f, _)| f)
}

pub fn parse_with(
    src: &str,
    origin: &str,
    opts: &ParseOptions,
) -> Result<(IndexFile, ParseReport)> {
    let mut file = IndexFile::new();

    for (lineno0, raw) in src.lines().enumerate() {
        let line = lineno0 + 1;
        // GROMACS strips comments before parsing, so a ';' inside a header is not protected.
        let content = match raw.find(';') {
            Some(i) => &raw[..i],
            None => raw,
        };
        let trimmed = content.trim();
        if trimmed.is_empty() {
            continue;
        }

        if trimmed.starts_with('[') {
            let inner = trimmed
                .strip_prefix('[')
                .and_then(|s| s.strip_suffix(']'))
                .ok_or_else(|| ParseError::UnterminatedHeader {
                    origin: origin.to_string(),
                    line,
                })?;
            let name = inner.trim();
            if name.is_empty() {
                return Err(ParseError::EmptyName {
                    origin: origin.to_string(),
                    line,
                }
                .into());
            }
            file.push(Group::new(name, Vec::new()));
            continue;
        }

        // An atom-index line. Column is 1-based within the raw line, for diagnostics.
        for (offset, token) in token_spans(content) {
            let col = offset + 1;
            let Some(group) = file.groups.last_mut() else {
                return Err(ParseError::AtomsBeforeGroup {
                    origin: origin.to_string(),
                    line,
                    col,
                }
                .into());
            };
            let atom: AtomId = match token.parse() {
                Ok(a) => a,
                Err(_) => {
                    // Distinguish "not a number at all" from "too big to fit".
                    let kind = if token.bytes().all(|b| b.is_ascii_digit()) {
                        ParseError::IndexTooLarge {
                            origin: origin.to_string(),
                            line,
                            col,
                            token: token.to_string(),
                        }
                    } else {
                        ParseError::BadToken {
                            origin: origin.to_string(),
                            line,
                            col,
                            token: token.to_string(),
                        }
                    };
                    return Err(kind.into());
                }
            };
            if atom == 0 {
                return Err(ParseError::ZeroIndex {
                    origin: origin.to_string(),
                    line,
                    col,
                }
                .into());
            }
            group.atoms.push(atom);
        }
    }

    let mut report = ParseReport::default();
    for g in &file.groups {
        if !g.is_sorted() {
            report.unsorted.push(g.name.clone());
        }
        if g.has_duplicates() {
            report.duplicated.push(g.name.clone());
        }
    }

    if opts.strict && !report.is_clean() {
        let mut msgs = Vec::new();
        if !report.unsorted.is_empty() {
            msgs.push(format!("unsorted group(s): {}", report.unsorted.join(", ")));
        }
        if !report.duplicated.is_empty() {
            msgs.push(format!(
                "group(s) with duplicate atoms: {}",
                report.duplicated.join(", ")
            ));
        }
        return Err(NdxError::Other(format!(
            "{origin}: {}\nhint: run `ndxed fmt --sort --dedup` to canonicalize",
            msgs.join("; ")
        )));
    }

    Ok((file, report))
}

/// Whitespace-delimited tokens with their byte offset in `line`.
fn token_spans(line: &str) -> impl Iterator<Item = (usize, &str)> {
    line.split_ascii_whitespace().map(move |tok| {
        // `tok` is a subslice of `line`, so pointer arithmetic gives the offset.
        let offset = tok.as_ptr() as usize - line.as_ptr() as usize;
        (offset, tok)
    })
}

pub fn parse_reader<R: Read>(mut r: R, origin: &str) -> Result<IndexFile> {
    let mut buf = String::new();
    r.read_to_string(&mut buf)
        .map_err(|e| NdxError::io(origin, e))?;
    parse_str(&buf, origin)
}

/// Read a `.ndx` file. `-` means stdin.
pub fn parse_path(p: &Path) -> Result<IndexFile> {
    parse_path_with(p, &ParseOptions::default()).map(|(f, _)| f)
}

pub fn parse_path_with(p: &Path, opts: &ParseOptions) -> Result<(IndexFile, ParseReport)> {
    let origin = p.display().to_string();
    let src = if p == Path::new("-") {
        let mut buf = String::new();
        std::io::stdin()
            .read_to_string(&mut buf)
            .map_err(|e| NdxError::io("<stdin>", e))?;
        buf
    } else {
        std::fs::read_to_string(p).map_err(|e| NdxError::io(p, e))?
    };
    let origin = if p == Path::new("-") {
        "<stdin>".to_string()
    } else {
        origin
    };
    parse_with(&src, &origin, opts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic() {
        let f = parse_str("[ System ]\n1 2 3\n[ SOL ]\n4 5\n", "t").unwrap();
        assert_eq!(f.len(), 2);
        assert_eq!(f.groups[0].name, "System");
        assert_eq!(f.groups[0].atoms, [1, 2, 3]);
        assert_eq!(f.groups[1].atoms, [4, 5]);
    }

    #[test]
    fn names_with_spaces_and_hyphens() {
        let f = parse_str("[ C-alpha ]\n1\n[ Water and ions ]\n2\n", "t").unwrap();
        assert_eq!(f.groups[0].name, "C-alpha");
        assert_eq!(f.groups[1].name, "Water and ions");
    }

    #[test]
    fn header_without_padding() {
        let f = parse_str("[System]\n1\n", "t").unwrap();
        assert_eq!(f.groups[0].name, "System");
    }

    #[test]
    fn comments_are_stripped() {
        let src = "; a leading comment\n[ A ]  ; trailing\n1 2 ; more\n; whole line\n3\n";
        let f = parse_str(src, "t").unwrap();
        assert_eq!(f.len(), 1);
        assert_eq!(f.groups[0].name, "A");
        assert_eq!(f.groups[0].atoms, [1, 2, 3]);
    }

    #[test]
    fn crlf_and_tabs() {
        let f = parse_str("[ A ]\r\n1\t2\r\n", "t").unwrap();
        assert_eq!(f.groups[0].atoms, [1, 2]);
    }

    #[test]
    fn no_trailing_newline() {
        let f = parse_str("[ A ]\n1 2", "t").unwrap();
        assert_eq!(f.groups[0].atoms, [1, 2]);
    }

    #[test]
    fn empty_group_is_legal() {
        let f = parse_str("[ A ]\n[ B ]\n1\n", "t").unwrap();
        assert_eq!(f.groups[0].atoms, [] as [AtomId; 0]);
        assert_eq!(f.groups[1].atoms, [1]);
    }

    #[test]
    fn duplicate_group_names_are_legal() {
        let f = parse_str("[ SOL ]\n1\n[ SOL ]\n2\n", "t").unwrap();
        assert_eq!(f.len(), 2);
        assert_eq!(f.find_by_name("SOL"), [0, 1]);
    }

    #[test]
    fn empty_file_is_legal() {
        assert_eq!(parse_str("", "t").unwrap().len(), 0);
    }

    #[test]
    fn report_flags_unsorted_and_dups() {
        let (_, r) = parse_with("[ A ]\n3 1 1\n[ B ]\n1 2\n", "t", &ParseOptions::default()).unwrap();
        assert_eq!(r.unsorted, ["A"]);
        assert_eq!(r.duplicated, ["A"]);
    }

    #[test]
    fn strict_rejects_unsorted() {
        let opts = ParseOptions { strict: true };
        assert!(parse_with("[ A ]\n3 1\n", "t", &opts).is_err());
    }

    #[test]
    fn err_atoms_before_group() {
        assert!(matches!(
            parse_str("1 2\n[ A ]\n", "t"),
            Err(NdxError::Parse(ParseError::AtomsBeforeGroup { line: 1, .. }))
        ));
    }

    #[test]
    fn err_zero_index() {
        assert!(matches!(
            parse_str("[ A ]\n1 0\n", "t"),
            Err(NdxError::Parse(ParseError::ZeroIndex { line: 2, col: 3, .. }))
        ));
    }

    #[test]
    fn err_bad_token() {
        match parse_str("[ A ]\n1 xy\n", "t") {
            Err(NdxError::Parse(ParseError::BadToken { line, col, token, .. })) => {
                assert_eq!((line, col, token.as_str()), (2, 3, "xy"));
            }
            other => panic!("expected BadToken, got {other:?}"),
        }
    }

    #[test]
    fn err_index_too_large() {
        assert!(matches!(
            parse_str("[ A ]\n99999999999999\n", "t"),
            Err(NdxError::Parse(ParseError::IndexTooLarge { .. }))
        ));
    }

    #[test]
    fn err_unterminated_header() {
        assert!(matches!(
            parse_str("[ A\n1\n", "t"),
            Err(NdxError::Parse(ParseError::UnterminatedHeader { line: 1, .. }))
        ));
    }

    #[test]
    fn err_empty_name() {
        assert!(matches!(
            parse_str("[  ]\n1\n", "t"),
            Err(NdxError::Parse(ParseError::EmptyName { line: 1, .. }))
        ));
    }
}
