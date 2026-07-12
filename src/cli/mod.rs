pub mod render;
pub mod run;

use std::path::PathBuf;

use clap::{ArgAction, Args, Parser, Subcommand, ValueEnum};

use ndx_editor::model::GroupRef;
use ndx_editor::ops::diff::Match;
use ndx_editor::ops::merge::OnConflict;
use ndx_editor::system::{BondMode, LoadOptions};

#[derive(Parser, Debug)]
#[command(
    name = "ndxed",
    version,
    about = "GROMACS .ndx index file editor",
    long_about = "Edit GROMACS index (.ndx) files from the shell, or interactively with a \
                  gmx make_ndx-style prompt.\n\n\
                  Mutating commands print the modified file to stdout by default and never touch \
                  the input file; pass -o FILE to write one.",
    propagate_version = true,
    args_conflicts_with_subcommands = true
)]
pub struct Cli {
    /// Suppress warnings on stderr.
    #[arg(short, long, global = true)]
    pub quiet: bool,

    /// Increase detail (repeatable).
    #[arg(short, long, global = true, action = ArgAction::Count)]
    pub verbose: u8,

    /// Open this file in the interactive editor. Same as `ndxed edit FILE`.
    pub file: Option<PathBuf>,

    #[command(subcommand)]
    pub cmd: Option<Cmd>,
}

/// Where a mutating command sends its output. The input file is never a target.
#[derive(Args, Debug, Clone)]
pub struct OutOpts {
    /// Write here instead of stdout.
    #[arg(short, long, value_name = "FILE")]
    pub out: Option<PathBuf>,

    /// Allow -o to point at the input file.
    #[arg(long)]
    pub force_overwrite: bool,

    /// Atom indices per line.
    #[arg(long, value_name = "N", default_value_t = 15)]
    pub per_line: usize,

    /// Minimum field width for an atom index.
    #[arg(long, value_name = "N", default_value_t = 4)]
    pub width: usize,
}

/// How `!` decides what "all atoms" means.
#[derive(Args, Debug, Clone)]
pub struct UniverseOpts {
    /// Total number of atoms in the system. Makes '!' exact.
    #[arg(long, value_name = "N", conflicts_with = "universe")]
    pub natoms: Option<u32>,

    /// Use this group as the universe for '!'.
    #[arg(long, value_name = "GROUP")]
    pub universe: Option<GroupRef>,
}

/// The structure and topology behind `name`, `element`, `bonded`, `within`, ...
#[derive(Args, Debug, Clone, Default)]
pub struct SystemOpts {
    /// Coordinates: a .gro file. Enables name/resname/resid/element/within.
    #[arg(short = 's', long, value_name = "FILE")]
    pub structure: Option<PathBuf>,

    /// Topology: a .top or .itp file. Enables real bonds, and type/molecule.
    #[arg(short = 'p', long, value_name = "FILE")]
    pub topology: Option<PathBuf>,

    /// Define a preprocessor symbol, as in an .mdp's `define = -DPOSRES`
    #[arg(short = 'D', value_name = "SYMBOL")]
    pub defines: Vec<String>,

    /// Search this directory for #include (repeatable; $GMXLIB is searched too)
    #[arg(short = 'I', value_name = "DIR")]
    pub include_dirs: Vec<PathBuf>,

    /// Where bonds come from
    #[arg(long, value_enum, default_value_t = BondMode::Auto)]
    pub bonds: BondMode,

    /// Apply the minimum-image convention in `within` (needs a box in the .gro)
    #[arg(long)]
    pub pbc: bool,
}

impl SystemOpts {
    pub fn load_options(&self) -> LoadOptions {
        LoadOptions {
            structure: self.structure.clone(),
            topology: self.topology.clone(),
            defines: self.defines.clone(),
            include_dirs: self.include_dirs.clone(),
            bonds: self.bonds,
            pbc: self.pbc,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum ListFormat {
    #[default]
    Table,
    Tsv,
    Json,
}

#[derive(Subcommand, Debug)]
pub enum Cmd {
    /// List the groups in an index file
    List {
        /// The .ndx file, or - for stdin
        file: PathBuf,

        /// Also show flags for unsorted groups and duplicate atoms
        #[arg(long)]
        long: bool,

        /// Show each group's atoms as compressed ranges
        #[arg(long)]
        ranges: bool,

        #[arg(long, value_enum, default_value_t = ListFormat::Table)]
        format: ListFormat,
    },

    /// Rename a group
    Rename {
        file: PathBuf,
        /// Group id (0-based) or name
        group: GroupRef,
        new_name: String,
        #[command(flatten)]
        out: OutOpts,
    },

    /// Add a group from a set expression, e.g. "0 & !1" or "Protein | SOL"
    #[command(long_about = "Add a group from a set expression.\n\n\
        Operators, tightest first:  !  (not)   &  (and)   \\  (difference)   |  (or)\n\
        Groups are named by 0-based id (`3`), by name (`SOL`), or by name and occurrence \
        (`SOL#2`). Quote a name that contains an operator: \"Protein & SOL\".\n\n\
        Examples:\n  \
          ndxed select sys.ndx \"0 & !1\"\n  \
          ndxed select sys.ndx \"Protein | SOL\" --name Solvated\n  \
          ndxed select sys.ndx \"(1 | 2) \\\\ 3\"")]
    Select {
        file: PathBuf,
        /// The expression
        expr: String,

        /// Name for the new group (default: a make_ndx-style name)
        #[arg(long)]
        name: Option<String>,

        /// Output only the new group
        #[arg(long)]
        keep_only: bool,

        /// Fail instead of creating an empty group
        #[arg(long)]
        error_on_empty: bool,

        #[command(flatten)]
        universe: UniverseOpts,
        #[command(flatten)]
        system: SystemOpts,
        #[command(flatten)]
        out: OutOpts,
    },

    /// Add a group from literal atom numbers, e.g. "1-10,15,20-30"
    Atoms {
        file: PathBuf,
        ranges: String,
        #[arg(long)]
        name: Option<String>,
        #[command(flatten)]
        out: OutOpts,
    },

    /// Delete groups
    Del {
        file: PathBuf,
        #[arg(required = true)]
        groups: Vec<GroupRef>,
        #[command(flatten)]
        out: OutOpts,
    },

    /// Keep only these groups, in the order given
    Keep {
        file: PathBuf,
        #[arg(required = true)]
        groups: Vec<GroupRef>,
        #[command(flatten)]
        out: OutOpts,
    },

    /// Add a group holding the first N atoms of a group (in file order)
    Head {
        file: PathBuf,
        group: GroupRef,
        n: usize,

        /// Take the last N instead
        #[arg(long)]
        tail: bool,

        #[arg(long)]
        name: Option<String>,

        /// Fail if the group has fewer than N atoms (default: take what there is)
        #[arg(long)]
        strict: bool,

        #[command(flatten)]
        out: OutOpts,
    },

    /// Split a group into several groups (atom order is preserved)
    Split {
        file: PathBuf,
        group: GroupRef,

        /// Split into K near-equal parts
        #[arg(long, value_name = "K", group = "how")]
        parts: Option<usize>,

        /// Split into chunks of K atoms
        #[arg(long, value_name = "K", group = "how")]
        size: Option<usize>,

        /// Cut after these 1-based positions within the group, e.g. 100,250
        #[arg(long, value_name = "POS", value_delimiter = ',', group = "how")]
        at: Option<Vec<usize>>,

        /// Split into one group per residue (-s) or per molecule (-p)
        #[arg(long, value_enum, group = "how")]
        by: Option<SplitBy>,

        /// Name prefix for the parts (default: the source group's name)
        #[arg(long)]
        prefix: Option<String>,

        /// Remove the source group
        #[arg(long)]
        replace: bool,

        #[command(flatten)]
        system: SystemOpts,
        #[command(flatten)]
        out: OutOpts,
    },

    /// Concatenate the groups of several index files
    Merge {
        /// Two or more .ndx files
        #[arg(required = true, num_args = 2..)]
        files: Vec<PathBuf>,

        /// What to do when two files define the same group name
        #[arg(long, value_enum, default_value_t = OnConflict::KeepBoth)]
        on_conflict: OnConflict,

        /// Prefix every group with its file's name, e.g. sysA:Protein
        #[arg(long)]
        prefix_file: bool,

        #[command(flatten)]
        out: OutOpts,
    },

    /// Compare two index files. Exits 1 if they differ, like diff(1)
    Diff {
        a: PathBuf,
        b: PathBuf,

        /// Also list the atoms that changed, as ranges
        #[arg(long)]
        atoms: bool,

        /// Pair groups up by name (default) or by position
        #[arg(long, value_enum, default_value_t = Match::Name)]
        by: Match,

        /// Always exit 0, even when the files differ
        #[arg(long)]
        no_exit_code: bool,
    },

    /// Rewrite an index file, optionally canonicalizing its groups
    Fmt {
        file: PathBuf,

        /// Sort each group's atoms
        #[arg(long)]
        sort: bool,

        /// Drop duplicate atoms within each group
        #[arg(long)]
        dedup: bool,

        #[command(flatten)]
        out: OutOpts,
    },

    /// Edit interactively, like gmx make_ndx. Type `help` at the prompt
    Edit {
        file: PathBuf,

        /// Where `q` saves. Without it, `q FILE` names the destination.
        #[arg(short, long, value_name = "FILE")]
        out: Option<PathBuf>,

        /// Allow saving over the input file
        #[arg(long)]
        force_overwrite: bool,

        /// Never write anything; just report what would happen
        #[arg(long)]
        dry_run: bool,

        /// Read plain lines instead of using the line editor
        #[arg(long)]
        no_readline: bool,

        #[command(flatten)]
        universe: UniverseOpts,
        #[command(flatten)]
        system: SystemOpts,
    },

    /// Summarize a structure/topology and check it against an index file
    Info {
        /// The .ndx file to check against (optional)
        file: Option<PathBuf>,
        #[command(flatten)]
        system: SystemOpts,
    },

    /// Print a shell completion script
    Completions { shell: clap_complete::Shell },
}

/// What `split --by` groups on.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum SplitBy {
    /// One group per residue id (needs -s conf.gro)
    Residue,
    /// One group per molecule instance (needs -p topol.top)
    Molecule,
}
