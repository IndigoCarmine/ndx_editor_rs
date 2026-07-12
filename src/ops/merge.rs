use std::path::Path;

use crate::model::{Group, IndexFile};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum OnConflict {
    /// Keep both. Duplicate names are legal in `.ndx`, and `#K` can pick between them.
    #[default]
    KeepBoth,
    /// Suffix the newcomer: `Protein` -> `Protein_2`.
    Rename,
    /// First one wins.
    Skip,
    /// Last one wins (the earlier group is replaced in place).
    Replace,
}

pub struct MergeReport {
    pub conflicts: Vec<String>,
}

/// Concatenate the groups of several index files into one, in file order.
pub fn merge<'a, I>(files: I, on_conflict: OnConflict, prefix_file: bool) -> (IndexFile, MergeReport)
where
    I: IntoIterator<Item = (&'a Path, &'a IndexFile)>,
{
    let mut out = IndexFile::new();
    let mut conflicts = Vec::new();

    for (path, f) in files {
        let prefix = prefix_file.then(|| stem(path));
        for g in &f.groups {
            let name = match &prefix {
                Some(p) => format!("{p}:{}", g.name),
                None => g.name.clone(),
            };
            let existing = out.find_by_name(&name);
            if existing.is_empty() {
                out.push(Group::new(name, g.atoms.clone()));
                continue;
            }

            conflicts.push(name.clone());
            match on_conflict {
                OnConflict::KeepBoth => {
                    out.push(Group::new(name, g.atoms.clone()));
                }
                OnConflict::Rename => {
                    let fresh = out.unique_name(&name);
                    out.push(Group::new(fresh, g.atoms.clone()));
                }
                OnConflict::Skip => {}
                OnConflict::Replace => {
                    let id = *existing.last().expect("non-empty");
                    out.groups[id] = Group::new(name, g.atoms.clone());
                }
            }
        }
    }

    (out, MergeReport { conflicts })
}

fn stem(p: &Path) -> String {
    p.file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| p.display().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn a() -> IndexFile {
        IndexFile {
            groups: vec![
                Group::new("System", vec![1, 2]),
                Group::new("Protein", vec![1]),
            ],
        }
    }

    fn b() -> IndexFile {
        IndexFile {
            groups: vec![
                Group::new("Protein", vec![9]),
                Group::new("Ligand", vec![3]),
            ],
        }
    }

    fn run(on_conflict: OnConflict, prefix: bool) -> IndexFile {
        let (fa, fb) = (a(), b());
        let files = [(Path::new("sysA.ndx"), &fa), (Path::new("sysB.ndx"), &fb)];
        merge(files, on_conflict, prefix).0
    }

    fn names(f: &IndexFile) -> Vec<String> {
        f.groups.iter().map(|g| g.name.clone()).collect()
    }

    #[test]
    fn keep_both_is_the_default() {
        let f = run(OnConflict::KeepBoth, false);
        assert_eq!(names(&f), ["System", "Protein", "Protein", "Ligand"]);
        assert_eq!(f.groups[1].atoms, [1]);
        assert_eq!(f.groups[2].atoms, [9]);
    }

    #[test]
    fn rename_suffixes_the_newcomer() {
        let f = run(OnConflict::Rename, false);
        assert_eq!(names(&f), ["System", "Protein", "Protein_2", "Ligand"]);
    }

    #[test]
    fn skip_keeps_the_first() {
        let f = run(OnConflict::Skip, false);
        assert_eq!(names(&f), ["System", "Protein", "Ligand"]);
        assert_eq!(f.groups[1].atoms, [1], "the original Protein survives");
    }

    #[test]
    fn replace_keeps_the_last() {
        let f = run(OnConflict::Replace, false);
        assert_eq!(names(&f), ["System", "Protein", "Ligand"]);
        assert_eq!(f.groups[1].atoms, [9], "the later Protein wins, in place");
    }

    #[test]
    fn prefix_file_avoids_conflicts_entirely() {
        let f = run(OnConflict::KeepBoth, true);
        assert_eq!(
            names(&f),
            ["sysA:System", "sysA:Protein", "sysB:Protein", "sysB:Ligand"]
        );
    }

    #[test]
    fn conflicts_are_reported() {
        let (fa, fb) = (a(), b());
        let files = [(Path::new("a.ndx"), &fa), (Path::new("b.ndx"), &fb)];
        let (_, report) = merge(files, OnConflict::Rename, false);
        assert_eq!(report.conflicts, ["Protein"]);
    }

    #[test]
    fn atom_order_is_preserved() {
        let fa = IndexFile {
            groups: vec![Group::new("Odd", vec![5, 1, 5])],
        };
        let files = [(Path::new("a.ndx"), &fa)];
        let (out, _) = merge(files, OnConflict::KeepBoth, false);
        assert_eq!(out.groups[0].atoms, [5, 1, 5]);
    }
}
