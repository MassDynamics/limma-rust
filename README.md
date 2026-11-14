# limma-rust

A pure-Rust port of the limma linear-model and empirical-Bayes core (`lmFit`, `contrasts.fit`,
`eBayes`, `topTable`, `decideTests`, `camera`) with R-faithful nmath and lowess. The crate is
`limma-core`.

It was split out of `md-limma` (`crates/limma-core` at md-limma `main` a4c3893) with its history.
The Python binding and the Mass Dynamics job layer stay in md-limma.

## Test

    cargo fmt --check && cargo clippy --all-targets -- -D warnings && cargo test

The tests compare against the R goldens in `tests-validation/` (inputs in `data/`, R's results in
`expected/`, the R scripts in `scripts/`; see its README). Set `MD_GOLDEN_CORPUS_DIR` to run them
against another corpus.
