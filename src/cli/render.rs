use std::io::Write;

use anstyle::{AnsiColor, Style};

use ndx_editor::atomset::{AtomSet, format_ranges};
use ndx_editor::model::IndexFile;
use ndx_editor::ops::diff::{Change, Diff};

use crate::cli::ListFormat;

const GREEN: Style = Style::new().fg_color(Some(anstyle::Color::Ansi(AnsiColor::Green)));
const RED: Style = Style::new().fg_color(Some(anstyle::Color::Ansi(AnsiColor::Red)));
const YELLOW: Style = Style::new().fg_color(Some(anstyle::Color::Ansi(AnsiColor::Yellow)));
const DIM: Style = Style::new().dimmed();
const BOLD: Style = Style::new().bold();

/// The group table, as `make_ndx` shows it on entry.
pub fn list(w: &mut impl Write, ndx: &IndexFile, long: bool, ranges: bool) -> std::io::Result<()> {
    if ndx.is_empty() {
        writeln!(w, "(no groups)")?;
        return Ok(());
    }
    let name_width = ndx
        .groups
        .iter()
        .map(|g| g.name.chars().count())
        .max()
        .unwrap_or(4)
        .clamp(4, 40);

    for (i, g) in ndx.groups.iter().enumerate() {
        write!(
            w,
            "{BOLD}{i:>3}{BOLD:#} {:<name_width$} {:>8} atoms",
            g.name,
            g.len()
        )?;
        if long {
            let mut flags = Vec::new();
            if !g.is_sorted() {
                flags.push("unsorted");
            }
            if g.has_duplicates() {
                flags.push("dups");
            }
            if !flags.is_empty() {
                write!(w, "   {YELLOW}{}{YELLOW:#}", flags.join(" "))?;
            }
        }
        writeln!(w)?;
        if ranges {
            let set = AtomSet::from_unsorted(g.atoms.clone());
            writeln!(w, "    {DIM}{}{DIM:#}", format_ranges(&set))?;
        }
    }
    Ok(())
}

pub fn list_tsv(w: &mut impl Write, ndx: &IndexFile) -> std::io::Result<()> {
    writeln!(w, "id\tname\tatoms")?;
    for (i, g) in ndx.groups.iter().enumerate() {
        writeln!(w, "{i}\t{}\t{}", g.name, g.len())?;
    }
    Ok(())
}

pub fn list_json(w: &mut impl Write, ndx: &IndexFile, ranges: bool) -> std::io::Result<()> {
    writeln!(w, "[")?;
    for (i, g) in ndx.groups.iter().enumerate() {
        let comma = if i + 1 == ndx.len() { "" } else { "," };
        write!(
            w,
            "  {{\"id\": {i}, \"name\": {}, \"atoms\": {}",
            json_string(&g.name),
            g.len()
        )?;
        if ranges {
            let set = AtomSet::from_unsorted(g.atoms.clone());
            write!(w, ", \"ranges\": {}", json_string(&format_ranges(&set)))?;
        }
        writeln!(w, "}}{comma}")?;
    }
    writeln!(w, "]")
}

fn json_string(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

pub fn list_any(
    w: &mut impl Write,
    ndx: &IndexFile,
    format: ListFormat,
    long: bool,
    ranges: bool,
) -> std::io::Result<()> {
    match format {
        ListFormat::Table => list(w, ndx, long, ranges),
        ListFormat::Tsv => list_tsv(w, ndx),
        ListFormat::Json => list_json(w, ndx, ranges),
    }
}

pub fn diff(
    w: &mut impl Write,
    d: &Diff,
    a_label: &str,
    b_label: &str,
    show_atoms: bool,
) -> std::io::Result<()> {
    writeln!(w, "{RED}--- {a_label}{RED:#}")?;
    writeln!(w, "{GREEN}+++ {b_label}{GREEN:#}")?;

    let mut same = 0usize;
    for e in &d.entries {
        match &e.change {
            Change::Same => same += 1,

            Change::Added => writeln!(
                w,
                "{GREEN}+ {} ({} atoms){GREEN:#}",
                e.name,
                e.len_b.unwrap_or(0)
            )?,

            Change::Removed => writeln!(
                w,
                "{RED}- {} ({} atoms){RED:#}",
                e.name,
                e.len_a.unwrap_or(0)
            )?,

            Change::Reordered => writeln!(
                w,
                "{YELLOW}~ {} (same atoms, different order){YELLOW:#}",
                e.name
            )?,

            Change::Changed { added, removed } => {
                writeln!(
                    w,
                    "{YELLOW}~ {} ({} -> {} atoms; +{} -{}){YELLOW:#}",
                    e.name,
                    e.len_a.unwrap_or(0),
                    e.len_b.unwrap_or(0),
                    added.len(),
                    removed.len()
                )?;
                if show_atoms {
                    if !removed.is_empty() {
                        writeln!(w, "    {RED}- {}{RED:#}", format_ranges(removed))?;
                    }
                    if !added.is_empty() {
                        writeln!(w, "    {GREEN}+ {}{GREEN:#}", format_ranges(added))?;
                    }
                }
            }
        }
    }

    if same > 0 {
        writeln!(w, "{DIM}= {same} group(s) unchanged{DIM:#}")?;
    }
    if d.is_empty() {
        writeln!(w, "{DIM}(no differences){DIM:#}")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndx_editor::model::Group;
    use ndx_editor::ops::diff::{Match, diff as make_diff};

    fn render<F>(f: F) -> String
    where
        F: FnOnce(&mut Vec<u8>) -> std::io::Result<()>,
    {
        let mut buf = Vec::new();
        f(&mut buf).unwrap();
        String::from_utf8(buf).unwrap()
    }

    fn ndx() -> IndexFile {
        IndexFile {
            groups: vec![
                Group::new("System", (1..=10).collect()),
                Group::new("Odd", vec![3, 1, 1]),
            ],
        }
    }

    #[test]
    fn table_shows_ids_names_and_counts() {
        let s = render(|w| list(w, &ndx(), false, false));
        assert!(s.contains("  0"));
        assert!(s.contains("System"));
        assert!(s.contains("10 atoms"));
    }

    #[test]
    fn long_flags_unsorted_and_dups() {
        let s = render(|w| list(w, &ndx(), true, false));
        assert!(s.contains("unsorted"));
        assert!(s.contains("dups"));
    }

    #[test]
    fn ranges_are_compressed() {
        let s = render(|w| list(w, &ndx(), false, true));
        assert!(s.contains("1-10"));
    }

    #[test]
    fn empty_file_says_so() {
        let s = render(|w| list(w, &IndexFile::new(), false, false));
        assert_eq!(s.trim(), "(no groups)");
    }

    #[test]
    fn tsv_is_machine_readable() {
        let s = render(|w| list_tsv(w, &ndx()));
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], "id\tname\tatoms");
        assert_eq!(lines[1], "0\tSystem\t10");
    }

    #[test]
    fn json_escapes_names() {
        let f = IndexFile {
            groups: vec![Group::new("we\"ird\\", vec![1])],
        };
        let s = render(|w| list_json(w, &f, false));
        assert!(s.contains(r#""name": "we\"ird\\""#), "{s}");
    }

    #[test]
    fn diff_marks_each_kind_of_change() {
        let a = IndexFile {
            groups: vec![
                Group::new("Same", vec![1]),
                Group::new("Gone", vec![2]),
                Group::new("Edited", vec![1, 2, 3]),
            ],
        };
        let b = IndexFile {
            groups: vec![
                Group::new("Same", vec![1]),
                Group::new("Edited", vec![2, 3, 4]),
                Group::new("Fresh", vec![9]),
            ],
        };
        let d = make_diff(&a, &b, Match::Name);
        let s = render(|w| diff(w, &d, "a.ndx", "b.ndx", true));
        assert!(s.contains("- Gone"));
        assert!(s.contains("+ Fresh"));
        assert!(s.contains("~ Edited"));
        assert!(s.contains("1 group(s) unchanged"));
        // --atoms itemizes the change.
        assert!(s.contains("- 1"));
        assert!(s.contains("+ 4"));
    }
}
