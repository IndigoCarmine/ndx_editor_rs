//! A GROMACS `.ndx` index-file editor.
//!
//! Two front ends sit on this library: a set of scriptable subcommands, and a `gmx make_ndx`-style
//! interactive REPL. Everything they can do lives here, so a future full-screen TUI can be added
//! without touching any logic.
//!
//! The central design decision is in [`model`]: a [`model::Group`] keeps its atoms in **file
//! order, with duplicates**, so reading and writing a file round-trips byte for byte. Set
//! arithmetic goes through [`atomset::AtomSet`], which is sorted and deduplicated — the same split
//! `make_ndx` makes.

pub mod atomset;
pub mod bonds;
pub mod error;
pub mod expr;
pub mod glob;
pub mod model;
pub mod ops;
pub mod parse;
pub mod spatial;
pub mod structure;
pub mod system;
pub mod topology;
pub mod universe;
pub mod write;

pub use error::{NdxError, Result};
pub use model::{AtomId, Group, GroupRef, IndexFile};
pub use parse::{parse_path, parse_str};
pub use write::{WriteOptions, to_string};
