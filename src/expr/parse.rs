use crate::error::ExprError;
use crate::expr::lex::{Keyword, Spanned, Tok, lex};
use crate::expr::{Expr, Pattern, Pred, Span};
use crate::model::GroupRef;

pub fn parse_expr(src: &str) -> Result<Expr, ExprError> {
    let toks = lex(src)?;
    if toks.is_empty() {
        return Err(ExprError::Empty);
    }
    let mut p = Parser {
        toks,
        pos: 0,
        end: src.len(),
    };
    let e = p.parse_or()?;
    if let Some(t) = p.peek_spanned() {
        return Err(ExprError::TrailingInput {
            span: t.span.clone(),
        });
    }
    Ok(e)
}

struct Parser {
    toks: Vec<Spanned>,
    pos: usize,
    end: usize,
}

impl Parser {
    fn peek_spanned(&self) -> Option<&Spanned> {
        self.toks.get(self.pos)
    }

    fn peek(&self) -> Option<&Tok> {
        self.toks.get(self.pos).map(|s| &s.tok)
    }

    fn peek_at(&self, n: usize) -> Option<&Tok> {
        self.toks.get(self.pos + n).map(|s| &s.tok)
    }

    fn bump(&mut self) -> Option<Spanned> {
        let t = self.toks.get(self.pos).cloned();
        if t.is_some() {
            self.pos += 1;
        }
        t
    }

    fn eat(&mut self, t: &Tok) -> bool {
        if self.peek() == Some(t) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    /// The span of the current token, or a zero-width span at end of input.
    fn here(&self) -> Span {
        match self.peek_spanned() {
            Some(s) => s.span.clone(),
            None => self.end..self.end,
        }
    }

    fn parse_or(&mut self) -> Result<Expr, ExprError> {
        let mut lhs = self.parse_and()?;
        while self.peek() == Some(&Tok::Pipe) {
            let op = self.bump().expect("peeked").span;
            let rhs = self.parse_and().map_err(|e| missing_operand(e, "|", op))?;
            lhs = Expr::Or(Box::new(lhs), Box::new(rhs));
        }
        Ok(lhs)
    }

    fn parse_and(&mut self) -> Result<Expr, ExprError> {
        let mut lhs = self.parse_unary()?;
        loop {
            let (is_diff, op) = match self.peek_spanned() {
                Some(s) if s.tok == Tok::Amp => (false, s.span.clone()),
                Some(s) if s.tok == Tok::Backslash => (true, s.span.clone()),
                _ => break,
            };
            self.pos += 1;
            let sym = if is_diff { "\\" } else { "&" };
            let rhs = self.parse_unary().map_err(|e| missing_operand(e, sym, op))?;
            lhs = if is_diff {
                Expr::Diff(Box::new(lhs), Box::new(rhs))
            } else {
                Expr::And(Box::new(lhs), Box::new(rhs))
            };
        }
        Ok(lhs)
    }

    fn parse_unary(&mut self) -> Result<Expr, ExprError> {
        match self.peek() {
            Some(Tok::Bang) => {
                let op = self.bump().expect("peeked").span;
                let inner = self.parse_unary().map_err(|e| missing_operand(e, "!", op))?;
                Ok(Expr::Not(Box::new(inner)))
            }

            Some(Tok::Kw(Keyword::Bonded)) => {
                let kw = self.bump().expect("peeked").span;
                // `bonded 2 of X` walks 2 bonds out; plain `bonded X` walks 1. Requiring `of` for
                // the depth form is what keeps `bonded 2` (the neighbours of group 2) unambiguous.
                let depth = match (self.peek(), self.peek_at(1)) {
                    (Some(Tok::Word { text, .. }), Some(Tok::Kw(Keyword::Of))) => {
                        match text.parse::<u8>() {
                            Ok(d) if d >= 1 => {
                                self.pos += 2;
                                d
                            }
                            _ => 1,
                        }
                    }
                    _ => 1,
                };
                let of = self
                    .parse_unary()
                    .map_err(|_| ExprError::MissingArgument {
                        keyword: "bonded",
                        span: kw.clone(),
                    })?;
                let span = kw.start..self.prev_end(kw.end);
                Ok(Expr::Bonded {
                    of: Box::new(of),
                    depth,
                    span,
                })
            }

            Some(Tok::Kw(Keyword::Within)) => {
                let kw = self.bump().expect("peeked").span;
                let radius = match self.bump() {
                    Some(Spanned {
                        tok: Tok::Word { text, .. },
                        span,
                    }) => text.trim().parse::<f32>().map_err(|_| ExprError::BadArgument {
                        msg: format!("`within` needs a radius in nm, got {text:?}"),
                        span,
                    })?,
                    _ => {
                        return Err(ExprError::MissingArgument {
                            keyword: "within",
                            span: kw,
                        });
                    }
                };
                if !self.eat(&Tok::Kw(Keyword::Of)) {
                    return Err(ExprError::ExpectedOf { span: self.here() });
                }
                let of = self
                    .parse_unary()
                    .map_err(|_| ExprError::MissingArgument {
                        keyword: "within",
                        span: kw.clone(),
                    })?;
                let span = kw.start..self.prev_end(kw.end);
                Ok(Expr::Within {
                    radius,
                    of: Box::new(of),
                    span,
                })
            }

            _ => self.parse_primary(),
        }
    }

    fn parse_primary(&mut self) -> Result<Expr, ExprError> {
        match self.peek() {
            Some(Tok::LParen) => {
                let open = self.bump().expect("peeked").span;
                let inner = self.parse_or()?;
                if !self.eat(&Tok::RParen) {
                    return Err(ExprError::UnbalancedParen { span: open });
                }
                Ok(inner)
            }

            Some(Tok::Kw(kw)) if is_predicate(*kw) => {
                let kw = *kw;
                let kw_span = self.bump().expect("peeked").span;
                let (arg, arg_span) = match self.bump() {
                    Some(Spanned {
                        tok: Tok::Word { text, .. },
                        span,
                    }) => (text, span),
                    _ => {
                        return Err(ExprError::MissingArgument {
                            keyword: keyword_name(kw),
                            span: kw_span,
                        });
                    }
                };
                let pred = build_pred(kw, &arg, &arg_span)?;
                Ok(Expr::Pred {
                    pred,
                    span: kw_span.start..arg_span.end,
                })
            }

            Some(Tok::Word { .. }) => {
                let Some(Spanned {
                    tok:
                        Tok::Word {
                            text,
                            occurrence,
                            quoted,
                        },
                    ..
                }) = self.bump()
                else {
                    unreachable!("just peeked a Word")
                };
                // Quoting turns off both readings: a bare `13` is group id 13 and a bare
                // `Fiber*` is a glob, but `"13"` and `"Fiber*"` are those literal names.
                Ok(Expr::Ref(match (quoted, text.parse::<usize>(), occurrence) {
                    (true, _, _) => GroupRef::Name {
                        name: text,
                        occurrence,
                    },
                    (false, Ok(id), None) => GroupRef::Id(id),
                    (false, _, None) => GroupRef::from_bare(text),
                    (false, _, Some(_)) => GroupRef::Name {
                        name: text,
                        occurrence,
                    },
                }))
            }

            _ => Err(ExprError::ExpectedGroupRef { span: self.here() }),
        }
    }

    /// End of the token we just consumed, for building a span that covers a whole sub-expression.
    fn prev_end(&self, fallback: usize) -> usize {
        self.toks
            .get(self.pos.wrapping_sub(1))
            .map(|s| s.span.end)
            .unwrap_or(fallback)
    }
}

/// A failure right after a binary operator reads better as "missing operand".
fn missing_operand(e: ExprError, op: &str, span: Span) -> ExprError {
    match e {
        ExprError::ExpectedGroupRef { .. } => ExprError::MissingOperand {
            op: op.to_string(),
            span,
        },
        other => other,
    }
}

fn is_predicate(kw: Keyword) -> bool {
    matches!(
        kw,
        Keyword::Name
            | Keyword::ResName
            | Keyword::ResId
            | Keyword::Element
            | Keyword::Chain
            | Keyword::Type
    )
}

fn keyword_name(kw: Keyword) -> &'static str {
    kw.as_str()
}

fn build_pred(kw: Keyword, arg: &str, span: &Span) -> Result<Pred, ExprError> {
    let words: Vec<&str> = arg.split_ascii_whitespace().collect();
    if words.is_empty() {
        return Err(ExprError::MissingArgument {
            keyword: keyword_name(kw),
            span: span.clone(),
        });
    }
    let pats = || words.iter().map(|w| Pattern::new(*w)).collect::<Vec<_>>();

    Ok(match kw {
        Keyword::Name => Pred::Name(pats()),
        Keyword::ResName => Pred::ResName(pats()),
        Keyword::Element => Pred::Element(pats()),
        Keyword::Type => Pred::Type(pats()),
        Keyword::Chain => {
            let mut cs = Vec::new();
            for w in &words {
                let mut it = w.chars();
                match (it.next(), it.next()) {
                    (Some(c), None) => cs.push(c),
                    _ => {
                        return Err(ExprError::BadArgument {
                            msg: format!("`chain` takes single characters, got {w:?}"),
                            span: span.clone(),
                        });
                    }
                }
            }
            Pred::Chain(cs)
        }
        Keyword::ResId => {
            let mut ranges = Vec::new();
            for tok in arg.split([',', ' ', '\t']).filter(|t| !t.is_empty()) {
                let bad = || ExprError::BadArgument {
                    msg: format!("`resid` takes numbers and ranges like `1-50`, got {tok:?}"),
                    span: span.clone(),
                };
                // Only split on a '-' that separates two numbers, so negative resids survive.
                let split = tok
                    .char_indices()
                    .skip(1)
                    .find(|(_, c)| *c == '-')
                    .map(|(i, _)| i);
                let (lo, hi) = match split {
                    Some(i) => {
                        let lo: i32 = tok[..i].parse().map_err(|_| bad())?;
                        let hi: i32 = tok[i + 1..].parse().map_err(|_| bad())?;
                        (lo, hi)
                    }
                    None => {
                        let n: i32 = tok.parse().map_err(|_| bad())?;
                        (n, n)
                    }
                };
                if lo > hi {
                    return Err(ExprError::BadArgument {
                        msg: format!("`resid` range {tok:?} runs backwards"),
                        span: span.clone(),
                    });
                }
                ranges.push((lo, hi));
            }
            if ranges.is_empty() {
                return Err(ExprError::MissingArgument {
                    keyword: "resid",
                    span: span.clone(),
                });
            }
            Pred::ResId(ranges)
        }
        Keyword::Bonded | Keyword::Within | Keyword::Of => unreachable!("not a predicate"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn id(n: usize) -> Expr {
        Expr::Ref(GroupRef::Id(n))
    }

    #[test]
    fn make_ndx_compatible() {
        assert_eq!(
            parse_expr("0 & !1").unwrap(),
            Expr::And(Box::new(id(0)), Box::new(Expr::Not(Box::new(id(1)))))
        );
    }

    #[test]
    fn and_binds_tighter_than_or() {
        // 0 | 1 & 2  ==  0 | (1 & 2)
        assert_eq!(
            parse_expr("0 | 1 & 2").unwrap(),
            Expr::Or(
                Box::new(id(0)),
                Box::new(Expr::And(Box::new(id(1)), Box::new(id(2))))
            )
        );
    }

    #[test]
    fn not_binds_tighter_than_and() {
        // !0 & 1  ==  (!0) & 1
        assert_eq!(
            parse_expr("!0 & 1").unwrap(),
            Expr::And(Box::new(Expr::Not(Box::new(id(0)))), Box::new(id(1)))
        );
    }

    #[test]
    fn parens_override() {
        assert_eq!(
            parse_expr("(0 | 1) & 2").unwrap(),
            Expr::And(
                Box::new(Expr::Or(Box::new(id(0)), Box::new(id(1)))),
                Box::new(id(2))
            )
        );
    }

    #[test]
    fn backslash_is_difference() {
        assert_eq!(
            parse_expr("0 \\ 1").unwrap(),
            Expr::Diff(Box::new(id(0)), Box::new(id(1)))
        );
    }

    #[test]
    fn or_is_left_associative() {
        assert_eq!(
            parse_expr("0 | 1 | 2").unwrap(),
            Expr::Or(
                Box::new(Expr::Or(Box::new(id(0)), Box::new(id(1)))),
                Box::new(id(2))
            )
        );
    }

    #[test]
    fn names_and_occurrences() {
        assert_eq!(parse_expr("Protein").unwrap(), Expr::Ref(GroupRef::name("Protein")));
        assert_eq!(
            parse_expr("SOL#2").unwrap(),
            Expr::Ref(GroupRef::Name {
                name: "SOL".into(),
                occurrence: Some(2)
            })
        );
        // A quoted all-digit name is a name, not an id.
        assert_eq!(parse_expr("\"13\"").unwrap(), Expr::Ref(GroupRef::name("13")));
    }

    #[test]
    fn a_bare_word_with_a_star_is_a_glob() {
        assert_eq!(
            parse_expr("Fiber*").unwrap(),
            Expr::Ref(GroupRef::Glob("Fiber*".into()))
        );
        assert_eq!(
            parse_expr("O?").unwrap(),
            Expr::Ref(GroupRef::Glob("O?".into()))
        );
        // Quoting turns it back into a literal name.
        assert_eq!(
            parse_expr("\"Fiber*\"").unwrap(),
            Expr::Ref(GroupRef::name("Fiber*"))
        );
        // No glob character, no glob.
        assert_eq!(
            parse_expr("C-alpha").unwrap(),
            Expr::Ref(GroupRef::name("C-alpha"))
        );
    }

    #[test]
    fn globs_compose_with_the_operators() {
        assert_eq!(
            parse_expr("Fiber* | fiber*").unwrap(),
            Expr::Or(
                Box::new(Expr::Ref(GroupRef::Glob("Fiber*".into()))),
                Box::new(Expr::Ref(GroupRef::Glob("fiber*".into())))
            )
        );
        assert!(parse_expr("Alkyl* & !O*").is_ok());
    }

    #[test]
    fn predicates() {
        assert_eq!(
            parse_expr("name CA").unwrap(),
            Expr::Pred {
                pred: Pred::Name(vec![Pattern::new("CA")]),
                span: 0..7
            }
        );
        assert!(matches!(
            parse_expr("name CA CB").unwrap(),
            Expr::Pred { pred: Pred::Name(ref p), .. } if p.len() == 2
        ));
        assert!(matches!(
            parse_expr("resid 1-50,60").unwrap(),
            Expr::Pred { pred: Pred::ResId(ref r), .. } if r == &[(1, 50), (60, 60)]
        ));
        assert!(matches!(
            parse_expr("chain A").unwrap(),
            Expr::Pred { pred: Pred::Chain(ref c), .. } if c == &['A']
        ));
    }

    #[test]
    fn resid_accepts_negative_numbers() {
        assert!(matches!(
            parse_expr("resid -3--1").unwrap(),
            Expr::Pred { pred: Pred::ResId(ref r), .. } if r == &[(-3, -1)]
        ));
    }

    #[test]
    fn the_target_expression() {
        // "the H atoms adjacent to Protein"
        let e = parse_expr("element H & bonded Protein").unwrap();
        match e {
            Expr::And(lhs, rhs) => {
                assert!(matches!(*lhs, Expr::Pred { pred: Pred::Element(_), .. }));
                assert!(matches!(*rhs, Expr::Bonded { depth: 1, .. }));
            }
            other => panic!("unexpected {other:?}"),
        }
    }

    #[test]
    fn bonded_depth_needs_of() {
        // `bonded 2` is the neighbours of group 2 ...
        assert!(matches!(
            parse_expr("bonded 2").unwrap(),
            Expr::Bonded { depth: 1, ref of, .. } if **of == id(2)
        ));
        // ... while `bonded 2 of X` walks two bonds out from X.
        assert!(matches!(
            parse_expr("bonded 2 of 1").unwrap(),
            Expr::Bonded { depth: 2, ref of, .. } if **of == id(1)
        ));
    }

    #[test]
    fn within() {
        assert!(matches!(
            parse_expr("within 0.5 of 1").unwrap(),
            Expr::Within { radius, .. } if (radius - 0.5).abs() < 1e-6
        ));
    }

    #[test]
    fn bonded_composes_with_parens_and_not() {
        assert!(parse_expr("bonded (1 | 2)").is_ok());
        assert!(parse_expr("!bonded 1").is_ok());
        assert!(parse_expr("bonded !1").is_ok());
    }

    #[test]
    fn errors_carry_spans() {
        assert!(matches!(parse_expr(""), Err(ExprError::Empty)));
        assert!(matches!(
            parse_expr("0 &"),
            Err(ExprError::MissingOperand { ref op, .. }) if op == "&"
        ));
        assert!(matches!(
            parse_expr("(0 | 1"),
            Err(ExprError::UnbalancedParen { .. })
        ));
        assert!(matches!(
            parse_expr("0 (1)"),
            // Two primaries in a row: the second is trailing input.
            Err(ExprError::TrailingInput { .. })
        ));
        assert!(matches!(
            parse_expr("within 0.5 (1)"),
            Err(ExprError::ExpectedOf { .. })
        ));
        assert!(matches!(
            parse_expr("within abc of 1"),
            Err(ExprError::BadArgument { .. })
        ));
        assert!(matches!(
            parse_expr("!"),
            Err(ExprError::MissingOperand { ref op, .. }) if op == "!"
        ));
    }

    #[test]
    fn trailing_input_span_points_at_the_extra_token() {
        match parse_expr("0 (1)") {
            Err(ExprError::TrailingInput { span }) => assert_eq!(span, 2..3),
            other => panic!("unexpected {other:?}"),
        }
    }

    /// Maximal-munch means an unquoted `0 1` is a group *named* "0 1", not two refs. That is the
    /// price of allowing `Water and ions`, and it is why numeric ids never contain spaces.
    #[test]
    fn a_bare_run_with_spaces_is_one_name() {
        assert_eq!(parse_expr("0 1").unwrap(), Expr::Ref(GroupRef::name("0 1")));
    }
}
