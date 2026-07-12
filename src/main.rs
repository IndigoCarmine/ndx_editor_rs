mod cli;
mod repl;

use clap::Parser;

fn main() {
    let args = cli::Cli::parse();
    std::process::exit(cli::run::run(args));
}
