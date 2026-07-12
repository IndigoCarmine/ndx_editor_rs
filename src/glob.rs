//! Shell-style glob matching, used for group names (`Fiber*`) and — once a structure file can be
//! loaded — for atom-name predicates (`name H*`).

/// A pattern with `*` (any run of characters) and `?` (exactly one).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Pattern(String);

impl Pattern {
    pub fn new(s: impl Into<String>) -> Self {
        Pattern(s.into())
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    pub fn matches(&self, s: &str) -> bool {
        matches(&self.0, s)
    }
}

/// Whether `s` contains anything a glob would treat specially.
pub fn is_pattern(s: &str) -> bool {
    s.contains(['*', '?'])
}

/// Iterative glob with backtracking on `*`.
pub fn matches(pat: &str, s: &str) -> bool {
    let (pat, s) = (pat.as_bytes(), s.as_bytes());
    let (mut p, mut i) = (0usize, 0usize);
    let (mut star, mut resume) = (usize::MAX, 0usize);

    while i < s.len() {
        if p < pat.len() && (pat[p] == b'?' || pat[p] == s[i]) {
            p += 1;
            i += 1;
        } else if p < pat.len() && pat[p] == b'*' {
            star = p;
            resume = i;
            p += 1;
        } else if star != usize::MAX {
            // Backtrack: let the last `*` swallow one more byte.
            p = star + 1;
            resume += 1;
            i = resume;
        } else {
            return false;
        }
    }

    while p < pat.len() && pat[p] == b'*' {
        p += 1;
    }
    p == pat.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn literals() {
        assert!(matches("CA", "CA"));
        assert!(!matches("CA", "CB"));
        assert!(!matches("CA", "CAA"));
    }

    #[test]
    fn star() {
        assert!(matches("Fiber*", "Fiber1"));
        assert!(matches("Fiber*", "Fiber"));
        assert!(!matches("Fiber*", "fiberA"), "matching is case-sensitive");
        assert!(matches("*", "anything"));
        assert!(matches("*", ""));
        assert!(matches("*A", "AlkylA"));
        assert!(matches("*yl*", "AlkylA"));
    }

    #[test]
    fn question_mark() {
        assert!(matches("?W", "OW"));
        assert!(!matches("?W", "HW1"));
        assert!(matches("Fiber?", "Fiber1"));
        assert!(!matches("Fiber?", "Fiber12"));
    }

    #[test]
    fn backtracking() {
        assert!(matches("C*A", "CBBA"));
        assert!(!matches("C*A", "CBBB"));
        assert!(matches("*a*b*", "xxaxxbxx"));
        assert!(!matches("*a*b*", "xxbxxaxx"));
    }

    #[test]
    fn empty_pattern_matches_only_empty() {
        assert!(matches("", ""));
        assert!(!matches("", "x"));
    }

    #[test]
    fn detecting_patterns() {
        assert!(is_pattern("O*"));
        assert!(is_pattern("Fiber?"));
        assert!(!is_pattern("C-alpha"));
        assert!(!is_pattern("Water and ions"));
    }
}
