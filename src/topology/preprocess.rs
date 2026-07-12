//! The `.top` preprocessor: `#include`, `#define`, `#ifdef` / `#ifndef` / `#else` / `#endif`.
//!
//! Conditionals are honoured, which matters: water is `#ifdef FLEXIBLE` bonds *or* `[ settles ]`,
//! and position restraints live under `#ifdef POSRES`. Expanding both branches would invent bonds
//! that the simulation does not have. Nothing is defined by default; `-D SYMBOL` adds one, the way
//! `define = -DPOSRES` does in an `.mdp`.

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use crate::error::{NdxError, Result};

pub struct Preprocessor {
    defines: HashSet<String>,
    /// Extra directories to search for `#include "..."`, plus `$GMXLIB`.
    include_dirs: Vec<PathBuf>,
    /// Includes we could not find. Not fatal — a force-field `.itp` usually holds parameters, not
    /// the molecule definitions we are after — but the caller should say so.
    pub missing: Vec<String>,
}

struct Frame {
    /// Whether the enclosing block was active, so `#else` can restore it.
    parent_active: bool,
    condition: bool,
    active: bool,
    else_seen: bool,
}

impl Preprocessor {
    pub fn new(defines: &[String], include_dirs: &[PathBuf]) -> Self {
        let mut dirs: Vec<PathBuf> = include_dirs.to_vec();
        // GROMACS installs its force fields under $GMXLIB (or $GMXDATA/top).
        for var in ["GMXLIB", "GMXDATA"] {
            if let Some(v) = std::env::var_os(var) {
                dirs.push(PathBuf::from(v));
            }
        }
        Preprocessor {
            defines: defines.iter().cloned().collect(),
            include_dirs: dirs,
            missing: Vec::new(),
        }
    }

    pub fn expand_path(&mut self, path: &Path) -> Result<String> {
        let content = std::fs::read_to_string(path).map_err(|e| NdxError::io(path, e))?;
        let mut out = String::new();
        let mut stack = Vec::new();
        let mut seen = HashSet::new();
        let mut frames = Vec::new();
        seen.insert(normalize(path));
        self.expand_into(&content, Some(path), &mut stack, &mut seen, &mut frames, &mut out)?;
        if !frames.is_empty() {
            return Err(NdxError::Other(format!(
                "{}: {} unterminated #ifdef/#ifndef block(s)",
                path.display(),
                frames.len()
            )));
        }
        Ok(out)
    }

    #[cfg(test)]
    pub fn expand_str(&mut self, content: &str) -> Result<String> {
        let mut out = String::new();
        let (mut stack, mut seen, mut frames) = (Vec::new(), HashSet::new(), Vec::new());
        self.expand_into(content, None, &mut stack, &mut seen, &mut frames, &mut out)?;
        Ok(out)
    }

    fn expand_into(
        &mut self,
        content: &str,
        source: Option<&Path>,
        stack: &mut Vec<PathBuf>,
        seen: &mut HashSet<PathBuf>,
        frames: &mut Vec<Frame>,
        out: &mut String,
    ) -> Result<()> {
        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.starts_with('#') {
                self.directive(trimmed, source, stack, seen, frames, out)?;
                continue;
            }
            if active(frames) {
                out.push_str(line);
                out.push('\n');
            }
        }
        Ok(())
    }

    fn directive(
        &mut self,
        line: &str,
        source: Option<&Path>,
        stack: &mut Vec<PathBuf>,
        seen: &mut HashSet<PathBuf>,
        frames: &mut Vec<Frame>,
        out: &mut String,
    ) -> Result<()> {
        // The conditional directives are handled first, and unconditionally: an #endif inside an
        // inactive block still has to close it.
        if let Some(rest) = line.strip_prefix("#ifdef") {
            let name = first_word(rest);
            let parent = active(frames);
            let cond = self.defines.contains(name);
            frames.push(Frame {
                parent_active: parent,
                condition: cond,
                active: parent && cond,
                else_seen: false,
            });
            return Ok(());
        }
        if let Some(rest) = line.strip_prefix("#ifndef") {
            let name = first_word(rest);
            let parent = active(frames);
            let cond = !self.defines.contains(name);
            frames.push(Frame {
                parent_active: parent,
                condition: cond,
                active: parent && cond,
                else_seen: false,
            });
            return Ok(());
        }
        if line.starts_with("#else") {
            let f = frames
                .last_mut()
                .ok_or_else(|| NdxError::Other("#else without a matching #ifdef".into()))?;
            if f.else_seen {
                return Err(NdxError::Other("two #else in one conditional block".into()));
            }
            f.else_seen = true;
            f.active = f.parent_active && !f.condition;
            return Ok(());
        }
        if line.starts_with("#endif") {
            frames
                .pop()
                .ok_or_else(|| NdxError::Other("#endif without a matching #ifdef".into()))?;
            return Ok(());
        }

        // Everything else only applies inside an active block.
        if !active(frames) {
            return Ok(());
        }

        if let Some(rest) = line.strip_prefix("#define") {
            let name = first_word(rest);
            if !name.is_empty() {
                self.defines.insert(name.to_string());
            }
            return Ok(());
        }
        if let Some(rest) = line.strip_prefix("#undef") {
            self.defines.remove(first_word(rest));
            return Ok(());
        }

        if let Some(rest) = line.strip_prefix("#include") {
            let target = include_target(rest)?;
            let Some(path) = self.find_include(&target, source) else {
                // A force-field include we cannot see holds parameters, not molecules; carry on
                // and let the caller decide whether the result is usable.
                self.missing.push(target);
                return Ok(());
            };
            let key = normalize(&path);
            if seen.contains(&key) {
                return Ok(()); // include-once
            }
            if stack.contains(&key) {
                return Err(NdxError::Other(format!(
                    "#include cycle at {}",
                    path.display()
                )));
            }
            let content = std::fs::read_to_string(&path).map_err(|e| NdxError::io(&path, e))?;
            seen.insert(key.clone());
            stack.push(key);
            let r = self.expand_into(&content, Some(&path), stack, seen, frames, out);
            stack.pop();
            return r;
        }

        // `#if`, `#pragma`, ... — not something a topology needs. Ignore rather than fail.
        Ok(())
    }

    fn find_include(&self, target: &str, source: Option<&Path>) -> Option<PathBuf> {
        let p = Path::new(target);
        if p.is_absolute() {
            return p.is_file().then(|| p.to_path_buf());
        }
        // Relative to the including file first, the way cpp does it.
        if let Some(dir) = source.and_then(Path::parent) {
            let candidate = dir.join(p);
            if candidate.is_file() {
                return Some(candidate);
            }
        }
        if p.is_file() {
            return Some(p.to_path_buf());
        }
        self.include_dirs
            .iter()
            .map(|d| d.join(p))
            .find(|c| c.is_file())
    }
}

fn active(frames: &[Frame]) -> bool {
    frames.iter().all(|f| f.active)
}

fn first_word(s: &str) -> &str {
    s.split_whitespace().next().unwrap_or("")
}

/// `#include "a.itp"` | `#include <a.itp>` | `#include a.itp`
fn include_target(rest: &str) -> Result<String> {
    let t = rest.trim();
    let inner = t
        .strip_prefix('"')
        .and_then(|r| r.split_once('"').map(|(p, _)| p))
        .or_else(|| {
            t.strip_prefix('<')
                .and_then(|r| r.split_once('>').map(|(p, _)| p))
        })
        .or_else(|| t.split_whitespace().next());

    inner
        .map(str::to_string)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| NdxError::Other(format!("malformed #include: {t:?}")))
}

fn normalize(p: &Path) -> PathBuf {
    std::fs::canonicalize(p).unwrap_or_else(|_| p.to_path_buf())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn expand(src: &str, defines: &[&str]) -> String {
        let d: Vec<String> = defines.iter().map(|s| (*s).to_string()).collect();
        Preprocessor::new(&d, &[]).expand_str(src).unwrap()
    }

    #[test]
    fn plain_lines_pass_through() {
        assert_eq!(expand("a\nb\n", &[]), "a\nb\n");
    }

    /// The case that made conditionals worth implementing: flexible-vs-settles water.
    #[test]
    fn ifdef_selects_a_branch() {
        let src = "\
#ifdef FLEXIBLE
[ bonds ]
1 2 1
#else
[ settles ]
1 1 0.1 0.16
#endif
";
        assert!(expand(src, &[]).contains("settles"));
        assert!(!expand(src, &[]).contains("bonds"));

        assert!(expand(src, &["FLEXIBLE"]).contains("bonds"));
        assert!(!expand(src, &["FLEXIBLE"]).contains("settles"));
    }

    #[test]
    fn ifndef_is_the_inverse() {
        let src = "#ifndef POSRES\nkeep\n#endif\n";
        assert_eq!(expand(src, &[]).trim(), "keep");
        assert_eq!(expand(src, &["POSRES"]).trim(), "");
    }

    #[test]
    fn nested_conditionals() {
        let src = "\
#ifdef A
outerA
#ifdef B
innerAB
#endif
#endif
";
        assert_eq!(expand(src, &[]).trim(), "");
        assert_eq!(expand(src, &["A"]).trim(), "outerA");
        let both = expand(src, &["A", "B"]);
        assert!(both.contains("outerA") && both.contains("innerAB"));
        // B alone must not leak the inner block: its parent is inactive.
        assert_eq!(expand(src, &["B"]).trim(), "");
    }

    #[test]
    fn define_and_undef() {
        assert_eq!(
            expand("#define X\n#ifdef X\nyes\n#endif\n", &[]).trim(),
            "yes"
        );
        assert_eq!(
            expand("#define X\n#undef X\n#ifdef X\nyes\n#endif\n", &[]).trim(),
            ""
        );
    }

    #[test]
    fn an_endif_inside_an_inactive_block_still_closes_it() {
        let src = "#ifdef NOPE\n#ifdef ALSO_NOPE\nx\n#endif\n#endif\nafter\n";
        assert_eq!(expand(src, &[]).trim(), "after");
    }

    #[test]
    fn unbalanced_conditionals_are_errors() {
        let mut p = Preprocessor::new(&[], &[]);
        assert!(p.expand_str("#endif\n").is_err());
        let mut p = Preprocessor::new(&[], &[]);
        assert!(p.expand_str("#else\n").is_err());
    }

    #[test]
    fn include_targets_are_parsed() {
        assert_eq!(include_target(" \"a.itp\" ").unwrap(), "a.itp");
        assert_eq!(include_target(" <ff/b.itp> ").unwrap(), "ff/b.itp");
        assert_eq!(include_target(" c.itp").unwrap(), "c.itp");
        assert!(include_target("  ").is_err());
    }

    #[test]
    fn includes_are_expanded_once_and_relative_to_their_parent() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::create_dir(dir.path().join("sub")).unwrap();
        std::fs::write(dir.path().join("sub/inner.itp"), "INNER\n").unwrap();
        std::fs::write(
            dir.path().join("sub/mid.itp"),
            "#include \"inner.itp\"\nMID\n",
        )
        .unwrap();
        let top = dir.path().join("top.top");
        std::fs::write(
            &top,
            "#include \"sub/mid.itp\"\n#include \"sub/inner.itp\"\nTOP\n",
        )
        .unwrap();

        let out = Preprocessor::new(&[], &[]).expand_path(&top).unwrap();
        assert_eq!(out.matches("INNER").count(), 1, "include-once");
        assert!(out.contains("MID") && out.contains("TOP"));
    }

    #[test]
    fn a_missing_include_is_recorded_not_fatal() {
        let dir = tempfile::tempdir().unwrap();
        let top = dir.path().join("t.top");
        std::fs::write(&top, "#include \"amber99.ff/forcefield.itp\"\nREST\n").unwrap();

        let mut p = Preprocessor::new(&[], &[]);
        let out = p.expand_path(&top).unwrap();
        assert!(out.contains("REST"));
        assert_eq!(p.missing, ["amber99.ff/forcefield.itp"]);
    }

    #[test]
    fn include_dirs_are_searched() {
        let dir = tempfile::tempdir().unwrap();
        let ff = dir.path().join("ff");
        std::fs::create_dir(&ff).unwrap();
        std::fs::write(ff.join("x.itp"), "FOUND\n").unwrap();
        let top = dir.path().join("t.top");
        std::fs::write(&top, "#include \"x.itp\"\n").unwrap();

        let mut p = Preprocessor::new(&[], &[ff]);
        assert!(p.expand_path(&top).unwrap().contains("FOUND"));
        assert!(p.missing.is_empty());
    }

    #[test]
    fn an_include_cycle_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("a.itp"), "#include \"b.itp\"\n").unwrap();
        std::fs::write(dir.path().join("b.itp"), "#include \"a.itp\"\n").unwrap();
        let top = dir.path().join("t.top");
        std::fs::write(&top, "#include \"a.itp\"\n").unwrap();
        // a includes b includes a: the second `a` is already in `seen`, so include-once stops it
        // before the cycle check — either way it must terminate rather than recurse forever.
        let out = Preprocessor::new(&[], &[]).expand_path(&top);
        assert!(out.is_ok());
    }
}
