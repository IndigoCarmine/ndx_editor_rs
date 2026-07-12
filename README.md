# ndxed

A GROMACS `.ndx` index-file editor. Two front ends over one library:

- **scriptable subcommands** — `ndxed list`, `ndxed select`, `ndxed diff`, …
- **an interactive editor** — `ndxed edit sys.ndx -o new.ndx`, modelled on `gmx make_ndx`

Everything works on the `.ndx` alone. Hand it a `.gro` and a `.top` as well and the expression
language grows atom names, residues, elements, force-field types, molecules and a real bond graph:

```sh
ndxed select sys.ndx -s conf.gro -p topol.top "element H & bonded Protein"
```

## Install (Ubuntu)

Needs **Rust 1.88 or newer**. Ubuntu's `apt install rustc` is too old (24.04 ships 1.75), so use
rustup — do *not* mix the two.

```sh
# 1. A C linker, which Rust needs to link the binary.
sudo apt update
sudo apt install -y build-essential curl

# 2. Rust, via rustup.
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"          # or just open a new shell

# 3. Build and install ndxed.
git clone <this-repo> ndx_editor_rs
cd ndx_editor_rs
cargo install --path .             # -> ~/.cargo/bin/ndxed

ndxed --version
```

`cargo install` puts `ndxed` in `~/.cargo/bin`, which rustup adds to your `PATH`. If `ndxed: command
not found`, either open a new shell or add it yourself:

```sh
echo 'export PATH="$HOME/.cargo/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc
```

**If you already have Rust** but it is too old, `rustup update` is enough. If Rust came from apt,
remove it first (`sudo apt remove rustc cargo`) so it does not shadow rustup's.

<details>
<summary>Other ways to install</summary>

```sh
# Just build, without installing: the binary lands in target/release/ndxed
cargo build --release

# System-wide, for every user
cargo build --release && sudo install -m755 target/release/ndxed /usr/local/bin/

# Without the interactive line editor (no rustyline; `ndxed edit` still works, reading plain lines)
cargo install --path . --no-default-features
```
</details>

### Shell completion

```sh
mkdir -p ~/.local/share/bash-completion/completions
ndxed completions bash > ~/.local/share/bash-completion/completions/ndxed
```

`zsh`, `fish`, `elvish` and `powershell` work too — pass the name instead of `bash`.

## Interactive

```
$ ndxed edit sys.ndx -o new.ndx

  0 System        22450 atoms
  1 Protein        1960 atoms
  2 SOL           20490 atoms

> 1 & !2
  3 Protein_&_!SOL : 1960 atoms
> name 3 Prot_dry
> head 1 100
  4 Protein_head100 : 100 atoms
> q
wrote 5 group(s) to new.ndx
```

`ndxed sys.ndx` (no subcommand) does the same thing.

**The input file is never overwritten.** `q` writes to `-o`, or to the path you give it (`q out.ndx`).
Ctrl-D behaves like `q`; `q!` quits without saving.

**Enter on an empty line re-prints the group table** — the thing you want most often after an edit.

Type `help` at the prompt for the full command list: `l`, `name`, `del`, `keep`, `a`, `head`,
`tail`, `split`, `merge`, `diff`, `undo`, `redo`, `w`, `q`, `q!`.

## Expressions

Operators, tightest binding first: `!` (not), `&` (and), `\` (difference), `|` (or). Parentheses work.

Groups are named by 0-based id (`3`), by name (`SOL`), by name and occurrence when a file has
duplicates (`SOL#2`), or by a **wildcard** (`Fiber*`, `O?`). Names may contain spaces and hyphens —
`C-alpha`, `Water and ions` — because `-` is deliberately *not* an operator and a bare name runs up
to the next operator. Quote a name that contains one: `"Protein & SOL"`.

```sh
ndxed select sys.ndx "0 & !1"                      # make_ndx-compatible
ndxed select sys.ndx "Protein | SOL" --name Solv
ndxed select sys.ndx "(1 | 2) \\ 3"
ndxed select sys.ndx "Protein & !SOL" --natoms 22450
```

## Wildcards

`*` matches any run of characters, `?` exactly one. **In an expression a wildcard means the union of
every group it matches**, so it composes with the operators:

```sh
ndxed select sys.ndx "Fiber*"              # Fiber1 ∪ Fiber2
ndxed select sys.ndx "Fiber* | fiber*"     # matching is case-sensitive
ndxed select sys.ndx "Alkyl* & !O*"
ndxed del    sys.ndx "O*"                  # delete OA and OB
ndxed keep   sys.ndx "*Monomer"
```

Quote your patterns — otherwise the *shell* tries to expand them against filenames first.

A pattern that matches no group is an error, not a silent no-op. Commands that can only act on one
group (`rename`, `split`, `--universe`) refuse a pattern that matches several, and say which:

```
$ ndxed rename sys.ndx "Fiber*" X
error: "Fiber*" matches 2 groups (Fiber1, Fiber2), but this takes exactly one
hint: name one of them, or use its id
```

Quoting *inside* an expression turns the wildcard off, so a group genuinely named `O*` is still
reachable — the same escape hatch that lets `"13"` mean the group named `13` rather than group id 13:

```sh
ndxed select sys.ndx 'O*'      # the pattern: every group starting with O
ndxed select sys.ndx '"O*"'    # the group literally called O*
```

`!` needs to know what "all atoms" means. With no structure file that is derived, in order:
`--natoms N` → `--universe GROUP` → the file's own `System` group → the union of all groups (with a
warning on stderr). `A & !B` is internally rewritten to `A \ B`, so the overwhelmingly common case
never needs a universe and never warns.

## Subcommands

| | |
|---|---|
| `ndxed list FILE [--long] [--ranges] [--format table\|tsv\|json]` | show the groups |
| `ndxed select FILE EXPR [--name N] [--keep-only]` | add a group from an expression |
| `ndxed atoms FILE "1-10,15,20-30"` | add a group from literal atom numbers |
| `ndxed rename FILE GROUP NEW_NAME` | rename |
| `ndxed del FILE GROUP...` / `keep FILE GROUP...` | delete / keep (keep also reorders; both take wildcards) |
| `ndxed head FILE GROUP N [--tail]` | first (or last) N atoms, **in file order** |
| `ndxed split FILE GROUP --parts K \| --size K \| --at 100,250 \| --by residue\|molecule` | split a group |
| `ndxed merge A.ndx B.ndx --on-conflict keep-both\|rename\|skip\|replace` | concatenate files |
| `ndxed diff A.ndx B.ndx [--atoms]` | compare; exits 1 if they differ, like `diff(1)` |
| `ndxed fmt FILE [--sort] [--dedup]` | canonicalize |
| `ndxed info [FILE] -s conf.gro -p topol.top` | summarize the system, check it against the .ndx |
| `ndxed completions SHELL` | shell completion script |

Mutating commands **print the modified file to stdout** and leave the input alone. `-o FILE` writes
a file instead (atomically — a crash never truncates it). `-` means stdin/stdout, so they compose:

```sh
ndxed select sys.ndx "0 & !1" | ndxed rename - 3 Prot_dry | ndxed list -
```

stdout is data; warnings and progress go to stderr. Exit codes: `0` ok, `1` `diff` found
differences, `2` bad arguments, `3` data error, `4` I/O error.

## Atom order

A group's atoms are stored **exactly as written — same order, duplicates kept** — because index
order is meaningful to some GROMACS tools. So `write(read(f)) == f` byte for byte, `head`/`tail`/
`split` mean "as written", and `diff` can tell you a group was *reordered*.

Set operations (`& | ! \`) are the exception: they always produce a sorted, deduplicated group,
exactly as `make_ndx` does. `list --long` flags any group that is unsorted or has duplicates, and
`ndxed fmt --sort --dedup` canonicalizes on request. Nothing is normalized behind your back.

## Coming from `gmx make_ndx`

| make_ndx | ndxed |
|---|---|
| `0 & !1` | same (in `edit`, or `ndxed select f.ndx "0 & !1"`) |
| `name 3 Foo` | same |
| `del 3` | same |
| `keep 3` | same |
| `l` | same |
| `q` (saves to `-o`) | same, but **refuses to overwrite the input** |
| `a 1-10` | same, plus `atomid 1-10` inside an expression |
| `splitres` | `split GROUP --by residue -s conf.gro` |
| `splitch` | `split GROUP --by molecule -p topol.top` (GROMACS has molecules, not chains) |
| `a CA` / `r 1-50` / `t opls_140` | `name CA` / `resid 1-50` / `type opls_140`, inside an expression |
| — | `head` / `tail`, `split --parts/--size/--at`, `merge`, `diff`, `undo`/`redo`, `\` for difference, parentheses, name references |

## Structure-dependent selection (`-s conf.gro`, `-p topol.top`)

Give `ndxed` a structure and/or a topology and the expression language grows per-atom predicates
and a bond graph. Everything below composes with the group algebra and the wildcards:

```sh
# The hydrogens attached to a group.
ndxed select sys.ndx -s conf.gro -p topol.top "element H & bonded Protein"

ndxed select sys.ndx -s conf.gro "resname SOL & within 0.5 of Protein"
ndxed select sys.ndx -s conf.gro "(name CA | name CB) & resid 1-50"
ndxed select sys.ndx -p topol.top "molecule SOL"
ndxed split  sys.ndx Protein --by residue  -s conf.gro     # make_ndx's splitres
ndxed split  sys.ndx SOL     --by molecule -p topol.top    # one group per water
ndxed info   sys.ndx -s conf.gro -p topol.top              # summary + consistency check
```

| predicate | needs | |
|---|---|---|
| `name CA`, `name H*` | `-s` | atom name (wildcards work; several values are OR'd) |
| `resname SOL` | `-s` | residue name |
| `resid 1-50,60` | `-s` | residue number |
| `element H` | `-s` | element (see below) |
| `type opls_140` | `-p` | force-field atom type |
| `molecule SOL` | `-p` | `[ moleculetype ]` name |
| `bonded X`, `bonded 2 of X` | bonds | atoms 1 (or N) bonds away from X |
| `within 0.5 of X` | `-s` | atoms within 0.5 nm of X |
| `atomid 116`, `atomid 1-10,15` | — | literal atom numbers; make_ndx's `a 1-10`, usable inside an expression |

`bonded X` and `within R of X` are *neighbourhoods*, so they contain X's own atoms when those are
bonded/close to each other — which is exactly what makes `element H & bonded Protein` pick up the
protein's own hydrogens. Write `bonded X & !X` for the strictly-outside version.

### Bonds

Real bonds come from the topology: `[ bonds ]`, and also **`[ constraints ]` and `[ settles ]`** —
without the latter a rigid water has no connectivity at all, since GROMACS never writes its O-H
bonds anywhere else.

With only `-s`, bonds are **estimated from interatomic distances** (covalent radii + 0.045 nm, via a
cell list). `ndxed` says so on stderr rather than pretending. `--bonds top|distance|auto` forces the
choice.

### Elements

Inferred from the atom name — **the first letter wins**, so `CA`, `CB`, `CG1` are all carbon and
`HW1`, `1HB` are hydrogen. A monatomic ion (`resname NA`, `atom NA`) keeps both letters. When a
topology is loaded, a **mass that matches a known element overrides the name**, which is what tells a
calcium ion apart from an alpha carbon. A mass that matches nothing — a coarse-grained bead, or
hydrogen-mass repartitioning — falls back to the name rather than inventing an element.

### The topology preprocessor

`#include` (searched relative to the including file, then `-I DIR`, then `$GMXLIB`), `#define`,
`#undef`, and **`#ifdef` / `#ifndef` / `#else` / `#endif`** are all honoured. That last part matters:
water is `#ifdef FLEXIBLE` bonds *or* `[ settles ]`, and expanding both branches would invent bonds
the simulation does not have. Nothing is defined by default; `-D SYMBOL` adds one, exactly like an
`.mdp`'s `define = -DPOSRES`.

A force-field `#include` that cannot be found is a **warning, not an error** — it usually holds only
parameters. But if a `[ molecules ]` line then names a molecule type nobody defined, that *is* an
error, and it tells you which include to go looking for.

### Without a system

The structure-aware forms still parse, and the failure names the flag that would have worked:

```
$ ndxed select sys.ndx "element H & bonded Protein"
error: `element` needs a structure file, which has not been loaded
hint: pass -s conf.gro
  element H & bonded Protein
  ^^^^^^^^^
```

`chain` is deliberately not supported: a `.gro` has no chain column, and GROMACS splits a system
into molecules rather than chains — so `molecule <name>` is what you want.

## Layout

`src/lib.rs` is the library (`ndx_editor`); `src/main.rs` is a ~10-line binary (`ndxed`). The
interactive editor is a thin front end over `repl::Session`, so a full-screen TUI could be added
beside it without touching any logic.

```
model.rs     Group / IndexFile / GroupRef      atomset.rs   sorted-unique set arithmetic
parse.rs     .ndx reader                       write.rs     .ndx writer (gmx byte format), atomic writes
glob.rs      wildcard matching                 universe.rs  what `!` complements against
expr/        lexer, parser, evaluator          ops/         select, split, merge, diff, …
structure/   .gro reader, element inference    topology/    .top reader + #include/#ifdef preprocessor
bonds.rs     bond graph (topology / distance)  spatial.rs   `within R of X`, with optional PBC
system.rs    loads and cross-checks the above  cli/         clap tree, rendering
repl/        session, commands, tty & piped drivers
```

## Tests

```sh
cargo test
cargo test --no-default-features   # builds without rustyline
cargo clippy --all-targets -- -D warnings
```

`tests/real_data.rs` runs against the real system in `test_data/` (512,328 atoms: a bundled fiber
in solvent, with velocities in the `.gro`, an `#ifdef INTER`-guarded
`[ intermolecular_interactions ]` block, and a `.ndx` covering only part of the box). That
directory is ~35 MB, so those tests **skip themselves when it is absent** rather than failing —
`cargo test` stays green on a bare clone.

The REPL's non-tty driver is a plain line loop, so it is tested without a pty:

```sh
printf 'l\n0 & !1\nname 3 Foo\nq out.ndx\n' | ndxed edit sys.ndx
```
