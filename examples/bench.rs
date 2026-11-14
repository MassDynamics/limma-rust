//! Times limma-rust on one dataset written by `tests-validation/scripts/bench_data.R` and
//! prints the same rows as `tests-validation/scripts/bench.R`, which it mirrors call for call:
//! impl, dataset, pipeline, stage, seconds (the median of `reps` runs).
//!
//! ```text
//! cargo run --release --example bench -- <dataset dir> <load|fit|block> <reps>
//! ```

use std::path::Path;
use std::time::Instant;
use std::{env, fs};

use limma_core::block::{duplicate_correlation, lm_fit_matrix_block};
use limma_core::contrasts::contrasts_fit;
use limma_core::ebayes::{ebayes, EBayesOptions};
use limma_core::fit::lm_fit_matrix;
use limma_core::toptable::{decide_tests_separate, top_table_t};

/// A `write.csv(..., row.names = FALSE)` table of numbers: `(nrow, ncol, column-major values)`,
/// with `NA` read as `NaN`.
fn read_csv(path: &Path) -> (usize, usize, Vec<f64>) {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut lines = text.lines();
    let ncol = lines.next().expect("header").split(',').count();
    let rows: Vec<Vec<f64>> = lines
        .map(|line| {
            line.split(',')
                .map(|v| {
                    if v == "NA" {
                        f64::NAN
                    } else {
                        v.parse().expect("number")
                    }
                })
                .collect()
        })
        .collect();
    let mut data = Vec::with_capacity(rows.len() * ncol);
    for j in 0..ncol {
        data.extend(rows.iter().map(|row| row[j]));
    }
    (rows.len(), ncol, data)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let dir = Path::new(&args[1]);
    let pipeline = args[2].as_str();
    let reps: usize = args[3].parse().expect("reps");
    let dataset = dir.file_name().unwrap().to_string_lossy().into_owned();

    let (ngenes, narrays, exprs) = read_csv(&dir.join("input_log2.csv"));
    let (_, ncoef, design) = read_csv(&dir.join("design.csv"));
    let (_, ncon, contrasts) = read_csv(&dir.join("contrasts.csv"));
    let block: Vec<usize> = read_csv(&dir.join("block.csv"))
        .2
        .iter()
        .map(|&b| b as usize - 1)
        .collect();

    let stage = |name: &str, f: &mut dyn FnMut()| {
        let mut secs: Vec<f64> = (0..reps)
            .map(|_| {
                let start = Instant::now();
                f();
                start.elapsed().as_secs_f64()
            })
            .collect();
        secs.sort_by(f64::total_cmp);
        println!(
            "Rust\t{dataset}\t{pipeline}\t{name}\t{:.6}",
            secs[secs.len() / 2]
        );
    };

    match pipeline {
        "fit" => {
            let mut fit = None;
            stage("lmFit", &mut || {
                fit = Some(lm_fit_matrix(&exprs, ngenes, narrays, &design, ncoef).unwrap())
            });
            let fit = fit.unwrap();
            let mut cfit = None;
            stage("contrasts.fit", &mut || {
                cfit = Some(contrasts_fit(&fit, &contrasts, ncon).unwrap())
            });
            let cfit = cfit.unwrap();
            let mut eb = None;
            stage("eBayes", &mut || {
                eb = Some(ebayes(&cfit, &EBayesOptions::default()).unwrap())
            });
            let eb = eb.unwrap();
            stage("topTable", &mut || {
                for j in 0..ncon {
                    top_table_t(&cfit, &eb, j, 0.95).unwrap();
                }
            });
            stage("decideTests", &mut || {
                decide_tests_separate(&cfit, &eb, 0.05);
            });
        }
        "block" => {
            let mut dc = None;
            stage("duplicateCorrelation", &mut || {
                dc = Some(
                    duplicate_correlation(&exprs, ngenes, narrays, &design, ncoef, &block, 0.15)
                        .unwrap(),
                )
            });
            let correlation = dc.unwrap().consensus_correlation;
            stage("lmFit(block)", &mut || {
                lm_fit_matrix_block(&exprs, ngenes, narrays, &design, ncoef, &block, correlation)
                    .unwrap();
            });
        }
        _ => {}
    }
}
