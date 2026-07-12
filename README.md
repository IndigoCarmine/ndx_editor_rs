# ndxed

A GROMACS `.ndx` index-file editor. Two front ends over one library:

- **scriptable subcommands** — `ndxed list`, `ndxed select`, `ndxed diff`, …
- **an interactive editor** — `ndxed edit sys.ndx -o new.ndx`, modelled on `gmx make_ndx`

No structure file is needed for anything below; `ndxed` works on the `.ndx` alone.

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
| `ndxed split FILE GROUP --parts K \| --size K \| --at 100,250` | split a group |
| `ndxed merge A.ndx B.ndx --on-conflict keep-both\|rename\|skip\|replace` | concatenate files |
| `ndxed diff A.ndx B.ndx [--atoms]` | compare; exits 1 if they differ, like `diff(1)` |
| `ndxed fmt FILE [--sort] [--dedup]` | canonicalize |
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
| `a 1-10` | same |
| `splitres` / `splitch` / `a CA` / `r 1-50` | not yet — they need a structure file (see below) |
| — | `head` / `tail`, `split --parts/--size/--at`, `merge`, `diff`, `undo`/`redo`, `\` for difference, parentheses, name references |

## Structure-dependent selection (not implemented yet)

The grammar already accepts the structure-aware forms, so that adding `.gro` / `.top` readers later
changes no syntax, no auto-generated names, and no existing expression:

```sh
ndxed select sys.ndx "element H & bonded Protein"   # the H atoms adjacent to Protein
ndxed select sys.ndx "resname SOL & within 0.5 of 1"
ndxed select sys.ndx "(name CA | name CB) & resid 1-50"
```

Today these parse and then fail with a specific message rather than a syntax error:

```
error: `bonded` needs bond information, which has not been loaded
hint: -p topol.top (real bonds) or -s conf.gro (distance estimate) is not supported yet
  element H & bonded Protein
              ^^^^^^^^^^^^^^
```

`src/system.rs` holds the `Structure` / `Topology` / `BondGraph` traits these will implement, and
`EvalCtx` / `Session` already take a `SystemCtx`. Filling them in means writing the readers and the
matching arms in `expr::eval` — nothing else moves.

## Layout

`src/lib.rs` is the library (`ndx_editor`); `src/main.rs` is a ~10-line binary (`ndxed`). The
interactive editor is a thin front end over `repl::Session`, so a full-screen TUI could be added
beside it without touching any logic.

```
model.rs     Group / IndexFile / GroupRef      atomset.rs  sorted-unique set arithmetic
parse.rs     .ndx reader                       write.rs    .ndx writer (gmx byte format), atomic writes
expr/        lexer, parser, evaluator          universe.rs what `!` complements against
ops/         select, split, merge, diff, …     system.rs   the .gro/.top extension point
cli/         clap tree, rendering              repl/        session, commands, tty & piped drivers
```

## Tests

```sh
cargo test
cargo test --no-default-features   # builds without rustyline
cargo clippy --all-targets -- -D warnings
```

The REPL's non-tty driver is a plain line loop, so it is tested without a pty:

```sh
printf 'l\n0 & !1\nname 3 Foo\nq out.ndx\n' | ndxed edit sys.ndx
```
