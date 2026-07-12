use crate::error::ExprError;
use crate::expr::Span;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Keyword {
    Name,
    ResName,
    ResId,
    Element,
    Chain,
    Type,
    Bonded,
    Within,
    Of,
}

impl Keyword {
    pub fn as_str(self) -> &'static str {
        match self {
            Keyword::Name => "name",
            Keyword::ResName => "resname",
            Keyword::ResId => "resid",
            Keyword::Element => "element",
            Keyword::Chain => "chain",
            Keyword::Type => "type",
            Keyword::Bonded => "bonded",
            Keyword::Within => "within",
            Keyword::Of => "of",
        }
    }
}

/// Longest first, so `resname` is not mistaken for `resid`'s prefix (and so on).
const KEYWORDS: &[(&str, Keyword)] = &[
    ("resname", Keyword::ResName),
    ("element", Keyword::Element),
    ("within", Keyword::Within),
    ("bonded", Keyword::Bonded),
    ("chain", Keyword::Chain),
    ("resnr", Keyword::ResId),
    ("resid", Keyword::ResId),
    ("name", Keyword::Name),
    ("type", Keyword::Type),
    ("of", Keyword::Of),
];

/// Characters that end a bare name.
const STOP: &[char] = &['(', ')', '|', '&', '!', '\\', '"', '#'];

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Tok {
    LParen,
    RParen,
    Amp,
    Pipe,
    Bang,
    Backslash,
    Kw(Keyword),
    /// A bare or quoted run of text: a group name, a number, or a predicate argument.
    Word {
        text: String,
        occurrence: Option<usize>,
        /// `"13"` is the *name* "13"; a bare `13` is group id 13.
        quoted: bool,
    },
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Spanned {
    pub tok: Tok,
    pub span: Span,
}

/// Tokenize an expression.
///
/// The interesting rule is the bare name: group names contain spaces and hyphens (`C-alpha`,
/// `Water and ions`), so a name cannot be a whitespace-delimited token. It is instead the longest
/// run up to the next operator, paren, quote — or the next *keyword at a word boundary*, which is
/// what lets `within 0.5 of Protein` stop the number before `of`.
pub fn lex(src: &str) -> Result<Vec<Spanned>, ExprError> {
    let b = src.as_bytes();
    let mut out = Vec::new();
    let mut i = 0usize;

    while i < b.len() {
        if b[i].is_ascii_whitespace() {
            i += 1;
            continue;
        }

        let start = i;
        let simple = match b[i] {
            b'(' => Some(Tok::LParen),
            b')' => Some(Tok::RParen),
            b'&' => Some(Tok::Amp),
            b'|' => Some(Tok::Pipe),
            b'!' => Some(Tok::Bang),
            b'\\' => Some(Tok::Backslash),
            _ => None,
        };
        if let Some(tok) = simple {
            i += 1;
            out.push(Spanned {
                tok,
                span: start..i,
            });
            continue;
        }

        if b[i] == b'"' {
            let Some(close) = src[i + 1..].find('"').map(|k| i + 1 + k) else {
                return Err(ExprError::UnterminatedQuote { span: i..src.len() });
            };
            let text = src[i + 1..close].to_string();
            i = close + 1;
            let occurrence = eat_occurrence(src, &mut i);
            out.push(Spanned {
                tok: Tok::Word {
                    text,
                    occurrence,
                    quoted: true,
                },
                span: start..i,
            });
            continue;
        }

        if let Some((kw, len)) = keyword_at(src, i) {
            i += len;
            out.push(Spanned {
                tok: Tok::Kw(kw),
                span: start..i,
            });
            continue;
        }

        // A bare word: run to the next stop char or the next keyword boundary.
        let end = bare_end(src, i);
        let text = src[i..end].trim_end().to_string();
        // The span must not cover the whitespace we just trimmed, or a caret would point past
        // the word it is blaming.
        let mut span_end = i + text.len();
        i = end;
        let occurrence = eat_occurrence(src, &mut i);
        if occurrence.is_some() {
            span_end = i;
        }
        out.push(Spanned {
            tok: Tok::Word {
                text,
                occurrence,
                quoted: false,
            },
            span: start..span_end,
        });
    }

    Ok(out)
}

/// A keyword starts at `i` only if it is followed by whitespace or `(` — so a group called
/// `named` lexes as a name, while `name CA` lexes as a predicate.
fn keyword_at(src: &str, i: usize) -> Option<(Keyword, usize)> {
    let rest = &src[i..];
    for (word, kw) in KEYWORDS {
        if let Some(after) = rest.strip_prefix(word) {
            match after.chars().next() {
                Some(c) if c.is_whitespace() || c == '(' => return Some((*kw, word.len())),
                _ => continue,
            }
        }
    }
    None
}

/// Where the bare word starting at `i` ends.
fn bare_end(src: &str, i: usize) -> usize {
    // A keyword cannot begin the word itself (`keyword_at` already ruled that out), so only look
    // for one after each whitespace run.
    let mut at_boundary = false;
    for (off, c) in src[i..].char_indices() {
        let j = i + off;
        if STOP.contains(&c) {
            return j;
        }
        if at_boundary && !c.is_whitespace() && keyword_at(src, j).is_some() {
            return j;
        }
        at_boundary = c.is_whitespace();
    }
    src.len()
}

/// `#2` immediately after a name picks the 2nd group with that name.
fn eat_occurrence(src: &str, i: &mut usize) -> Option<usize> {
    let rest = &src[*i..];
    let digits = rest.strip_prefix('#')?;
    let n: String = digits.chars().take_while(char::is_ascii_digit).collect();
    let k: usize = n.parse().ok()?;
    if k == 0 {
        return None;
    }
    *i += 1 + n.len();
    Some(k)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toks(src: &str) -> Vec<Tok> {
        lex(src).unwrap().into_iter().map(|s| s.tok).collect()
    }

    fn word(text: &str) -> Tok {
        Tok::Word {
            text: text.into(),
            occurrence: None,
            quoted: false,
        }
    }

    fn quoted(text: &str) -> Tok {
        Tok::Word {
            text: text.into(),
            occurrence: None,
            quoted: true,
        }
    }

    #[test]
    fn make_ndx_expression() {
        assert_eq!(toks("0 & !1"), [word("0"), Tok::Amp, Tok::Bang, word("1")]);
    }

    #[test]
    fn names_with_spaces_and_hyphens() {
        assert_eq!(
            toks("Water and ions | C-alpha"),
            [word("Water and ions"), Tok::Pipe, word("C-alpha")]
        );
    }

    #[test]
    fn hyphen_is_not_an_operator() {
        assert_eq!(toks("C-alpha"), [word("C-alpha")]);
    }

    #[test]
    fn backslash_is_difference() {
        assert_eq!(toks("A \\ B"), [word("A"), Tok::Backslash, word("B")]);
    }

    #[test]
    fn quoted_name_may_contain_operators() {
        assert_eq!(toks("\"Protein & SOL\""), [quoted("Protein & SOL")]);
    }

    /// A quoted all-digit token is a name, not an id. Without the `quoted` flag the parser
    /// could not tell them apart.
    #[test]
    fn quoting_survives_lexing() {
        assert_eq!(toks("\"13\""), [quoted("13")]);
        assert_eq!(toks("13"), [word("13")]);
    }

    #[test]
    fn occurrence_suffix() {
        assert_eq!(
            toks("SOL#2"),
            [Tok::Word {
                text: "SOL".into(),
                occurrence: Some(2),
                quoted: false
            }]
        );
    }

    #[test]
    fn keywords_need_a_following_space() {
        assert_eq!(toks("name CA"), [Tok::Kw(Keyword::Name), word("CA")]);
        // A group that merely starts with a keyword is still a name.
        assert_eq!(toks("named"), [word("named")]);
        assert_eq!(toks("nameless | 1"), [word("nameless"), Tok::Pipe, word("1")]);
    }

    #[test]
    fn a_keyword_boundary_ends_a_bare_word() {
        // This is what makes `within 0.5 of X` work: "0.5" must not swallow "of X".
        assert_eq!(
            toks("within 0.5 of Protein"),
            [
                Tok::Kw(Keyword::Within),
                word("0.5"),
                Tok::Kw(Keyword::Of),
                word("Protein")
            ]
        );
    }

    #[test]
    fn bonded_takes_an_expression() {
        assert_eq!(
            toks("element H & bonded Protein"),
            [
                Tok::Kw(Keyword::Element),
                word("H"),
                Tok::Amp,
                Tok::Kw(Keyword::Bonded),
                word("Protein")
            ]
        );
    }

    #[test]
    fn keyword_before_paren() {
        assert_eq!(
            toks("bonded (1 | 2)"),
            [
                Tok::Kw(Keyword::Bonded),
                Tok::LParen,
                word("1"),
                Tok::Pipe,
                word("2"),
                Tok::RParen
            ]
        );
    }

    #[test]
    fn resname_is_not_read_as_resid() {
        assert_eq!(toks("resname SOL"), [Tok::Kw(Keyword::ResName), word("SOL")]);
    }

    #[test]
    fn unterminated_quote() {
        assert!(matches!(
            lex("\"abc"),
            Err(ExprError::UnterminatedQuote { .. })
        ));
    }

    #[test]
    fn spans_point_at_the_source() {
        let ts = lex("0 & !1").unwrap();
        assert_eq!(ts[0].span, 0..1, "the word span must exclude the trailing space");
        assert_eq!(ts[1].span, 2..3);
        assert_eq!(ts[3].span, 5..6);
    }
}
