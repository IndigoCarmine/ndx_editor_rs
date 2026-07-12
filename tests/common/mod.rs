#![allow(dead_code)] // shared across test binaries; each uses a subset

use std::path::{Path, PathBuf};

use assert_cmd::Command;

pub fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

pub fn ndxed() -> Command {
    Command::cargo_bin("ndxed").expect("the ndxed binary is built")
}

/// Read a fixture's bytes, so a test can prove a command left the input alone.
pub fn bytes(p: &Path) -> Vec<u8> {
    std::fs::read(p).expect("fixture is readable")
}

/// Parse whatever a command printed to stdout back into an index file.
///
/// Asserting on the parsed structure rather than on bytes is the right default: it keeps the
/// tests about behaviour, and lets the writer's formatting be pinned in exactly one place.
pub fn parse_out(stdout: &[u8]) -> ndx_editor::IndexFile {
    let s = std::str::from_utf8(stdout).expect("output is UTF-8");
    ndx_editor::parse_str(s, "<stdout>").expect("output is a valid .ndx")
}

pub fn names(f: &ndx_editor::IndexFile) -> Vec<String> {
    f.groups.iter().map(|g| g.name.clone()).collect()
}
