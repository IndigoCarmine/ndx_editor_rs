use std::collections::VecDeque;
use std::path::{Path, PathBuf};

use ndx_editor::error::{NdxError, Result};
use ndx_editor::model::IndexFile;
use ndx_editor::system::SystemCtx;
use ndx_editor::universe::UniverseSpec;
use ndx_editor::write::{Sink, WriteOptions, write_atomic};

const UNDO_DEPTH: usize = 20;

pub struct Session {
    pub ndx: IndexFile,
    /// The file we opened. Never a write target unless `force_overwrite`.
    pub source: PathBuf,
    pub out: Option<PathBuf>,
    pub universe: UniverseSpec,
    /// Empty for now: `.gro` / `.top` loading is the next milestone. See `ndx_editor::system`.
    pub system: SystemCtx,
    pub write_opts: WriteOptions,
    pub dirty: bool,
    pub dry_run: bool,
    pub force_overwrite: bool,

    undo: VecDeque<IndexFile>,
    redo: Vec<IndexFile>,
}

impl Session {
    pub fn open(
        source: PathBuf,
        ndx: IndexFile,
        out: Option<PathBuf>,
        universe: UniverseSpec,
        dry_run: bool,
        force_overwrite: bool,
    ) -> Self {
        Session {
            ndx,
            source,
            out,
            universe,
            system: SystemCtx::default(),
            write_opts: WriteOptions::default(),
            dirty: false,
            dry_run,
            force_overwrite,
            undo: VecDeque::new(),
            redo: Vec::new(),
        }
    }

    /// Snapshot the whole index file before a mutation.
    ///
    /// A whole `.ndx` is at most a few MB and the stack is bounded, so this is cheaper than an
    /// inverse-operation log is to get right — `del`, which renumbers every id, is trivially
    /// correct this way.
    pub fn snapshot(&mut self) {
        self.undo.push_back(self.ndx.clone());
        if self.undo.len() > UNDO_DEPTH {
            self.undo.pop_front();
        }
        self.redo.clear();
    }

    pub fn mark_dirty(&mut self) {
        self.dirty = true;
    }

    /// Undo the snapshot a failed command took, without making it look like an edit.
    ///
    /// A command that errors must leave no trace: otherwise `undo` would first have to unwind the
    /// failure before it reached the user's last real edit.
    pub fn rollback(&mut self) {
        if let Some(prev) = self.undo.pop_back() {
            self.ndx = prev;
        }
    }

    pub fn undo(&mut self) -> Result<()> {
        let prev = self
            .undo
            .pop_back()
            .ok_or_else(|| NdxError::Other("nothing to undo".into()))?;
        self.redo.push(std::mem::replace(&mut self.ndx, prev));
        self.dirty = true;
        Ok(())
    }

    pub fn redo(&mut self) -> Result<()> {
        let next = self
            .redo
            .pop()
            .ok_or_else(|| NdxError::Other("nothing to redo".into()))?;
        self.undo.push_back(std::mem::replace(&mut self.ndx, next));
        self.dirty = true;
        Ok(())
    }

    /// Write to `to`, or to the session's output path. Never writes over the input file.
    pub fn save(&mut self, to: Option<&Path>) -> Result<PathBuf> {
        let target = to
            .map(Path::to_path_buf)
            .or_else(|| self.out.clone())
            .ok_or(NdxError::NoOutputPath)?;

        if !self.force_overwrite && same_file(&target, &self.source) {
            return Err(NdxError::WouldOverwriteInput { path: target });
        }

        if !self.dry_run {
            if target == Path::new("-") {
                Sink::Stdout.commit(&self.ndx, &self.write_opts, &[], true)?;
            } else {
                write_atomic(&target, &self.ndx, &self.write_opts)?;
            }
        }

        // Remember it, so a later bare `w` / `q` goes to the same place.
        self.out = Some(target.clone());
        self.dirty = false;
        Ok(target)
    }
}

fn same_file(a: &Path, b: &Path) -> bool {
    match (a.canonicalize(), b.canonicalize()) {
        (Ok(a), Ok(b)) => a == b,
        _ => a == b,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndx_editor::model::Group;

    fn session() -> Session {
        let ndx = IndexFile {
            groups: vec![Group::new("A", vec![1, 2])],
        };
        Session::open(
            PathBuf::from("in.ndx"),
            ndx,
            None,
            UniverseSpec::Auto,
            false,
            false,
        )
    }

    #[test]
    fn undo_restores_the_previous_state() {
        let mut s = session();
        s.snapshot();
        s.ndx.push(Group::new("B", vec![3]));
        assert_eq!(s.ndx.len(), 2);
        s.undo().unwrap();
        assert_eq!(s.ndx.len(), 1);
    }

    #[test]
    fn redo_reapplies() {
        let mut s = session();
        s.snapshot();
        s.ndx.push(Group::new("B", vec![3]));
        s.undo().unwrap();
        s.redo().unwrap();
        assert_eq!(s.ndx.len(), 2);
    }

    #[test]
    fn a_new_edit_clears_the_redo_stack() {
        let mut s = session();
        s.snapshot();
        s.ndx.push(Group::new("B", vec![3]));
        s.undo().unwrap();
        s.snapshot();
        s.ndx.push(Group::new("C", vec![4]));
        assert!(s.redo().is_err());
    }

    #[test]
    fn undo_on_an_empty_stack_errors() {
        assert!(session().undo().is_err());
    }

    #[test]
    fn undo_survives_a_delete_that_renumbers() {
        let mut s = session();
        s.ndx.push(Group::new("B", vec![3]));
        s.snapshot();
        s.ndx.remove_many(&[0]);
        assert_eq!(s.ndx.groups[0].name, "B");
        s.undo().unwrap();
        assert_eq!(s.ndx.groups[0].name, "A");
        assert_eq!(s.ndx.groups[1].name, "B");
    }

    #[test]
    fn the_undo_stack_is_bounded() {
        let mut s = session();
        for _ in 0..(UNDO_DEPTH + 5) {
            s.snapshot();
            s.ndx.push(Group::new("X", vec![9]));
        }
        assert_eq!(s.undo.len(), UNDO_DEPTH);
    }

    #[test]
    fn save_without_a_target_errors() {
        let mut s = session();
        assert!(matches!(s.save(None), Err(NdxError::NoOutputPath)));
    }

    #[test]
    fn save_refuses_to_overwrite_the_input() {
        let mut s = session();
        let same = s.source.clone();
        assert!(matches!(
            s.save(Some(&same)),
            Err(NdxError::WouldOverwriteInput { .. })
        ));
    }

    #[test]
    fn saving_sets_the_default_target() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.ndx");
        let mut s = session();
        s.dirty = true;
        s.save(Some(&out)).unwrap();
        assert_eq!(s.out.as_deref(), Some(out.as_path()));
        assert!(!s.dirty);
        // A later bare save reuses it.
        s.dirty = true;
        assert_eq!(s.save(None).unwrap(), out);
    }

    #[test]
    fn dry_run_writes_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.ndx");
        let mut s = session();
        s.dry_run = true;
        s.save(Some(&out)).unwrap();
        assert!(!out.exists());
    }
}
