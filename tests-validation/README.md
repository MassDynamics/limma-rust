# tests-validation

Parity of limma-rust with R limma (3.68.5, R 4.5.0, statmod 1.5.2), laid out as the team's
`tests-validation` convention: inputs, the R that produced the reference, the reference results,
the Rust results of the current run, and a report.

| Folder | Holds |
|---|---|
| `data/` | the inputs each case is computed from: `input_log2.csv`, `design.csv`, `contrasts.csv`, `block.csv`, `case.json` |
| `scripts/` | the R that wrote `expected/` (and `data/` for the block cases), and `report`, which builds the parity report |
| `expected/` | R's results: `scalar/` (nmath and lowess grids, inputs in the same tables), `matrix/<case>/` (lmFit, contrasts.fit, eBayes, topTable, decideTests), `block/<case>/` (duplicateCorrelation and the block path) |
| `out/` | the current run, gitignored: one row per compared file in `comparisons/` (written by every `cargo test`), and `parity.tsv`, `parity.md` and `cargo-test.log` from `scripts/report` |
| `report/` | the committed parity report, `parity.tsv` and `parity.md`: per compared file, the values checked, the largest relative error and where, and how many missed tolerance |

## Run

    cargo test                            # every golden is a test; fails past tolerance
    tests-validation/scripts/report       # the parity report, in Docker
    tests-validation/scripts/report --update

Every golden is a test, so a change that moves any result past its tolerance fails `cargo test`
and CI. Tolerances are in `tests/common/mod.rs` and the individual `tests/*.rs`. Set
`MD_GOLDEN_CORPUS_DIR` to run the scalar and matrix tests against the full corpus instead.

`scripts/report` runs the same tests in `rust:1.95.0` (linux/amd64, release) and compares the
per-file numbers with the committed `report/parity.tsv`. It fails when a value misses tolerance,
and also when any number moves at all, even inside tolerance. Review the diff it prints; if the
change is intended, rerun with `--update` and commit `report/` with the change. The numbers are
pinned to the image because the last digits differ between macOS arm64 and linux amd64 (debug and
release agree on amd64).

## Regenerate

- `scripts/block_goldens.R`, from the repo root in the corpus R image:

      docker run --rm --platform linux/amd64 -e LANG=C.UTF-8 -e LC_COLLATE=C -v $PWD:/w -w /w \
        md-flexi-r45-limma368 Rscript tests-validation/scripts/block_goldens.R

- `scripts/scalar_goldens.R` is a copy of `MDFlexiComparisons/data-raw/golden-corpus/scalar_goldens.R`.
  It reads the bojkova2020 intensities from the full MDFlexiComparisons golden corpus, so it does
  not run from this package alone. It writes that corpus's mixed layout (inputs and results in one
  case directory), which then has to be split into `data/` and `expected/`.

## Data

All inputs are small and committed (30 MB in total). The matrix cases are a 600-protein slice of
the public bojkova2020 dataset; the block cases are synthetic or built from that slice.
