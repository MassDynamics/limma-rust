//! The block path against `tests/golden/block/<case>/` (written by `block_goldens.R` there):
//! `duplicateCorrelation(block = )` and `lmFit(block = , correlation = )`. The fit is fed R's
//! consensus correlation so it is checked on its own. Tolerances are the matrix goldens' (logFC
//! abs 1e-10, everything else rel 1e-8); the correlations are rel 1e-8 too.

mod common;

use std::path::{Path, PathBuf};

use common::matrix::Matrix;
use common::{Report, Tol};
use limma_core::block::{duplicate_correlation, lm_fit_matrix_block};

const ABS_1E10: Tol = Tol {
    rel: 0.0,
    abs_floor: 1e-10,
};
const REL_1E8: Tol = Tol {
    rel: 1e-8,
    abs_floor: 0.0,
};

const CASES: [&str; 4] = [
    "techrep_synth",
    "techrep_synth_na",
    "bojkova_pairs",
    "techrep_mixed",
];

fn case_dir(case: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/block")
        .join(case)
}

/// Block labels coded in sorted order, as `factor()` levels them (the goldens run in the C
/// collation).
fn read_block(dir: &Path) -> Vec<usize> {
    let mut reader = csv::Reader::from_path(dir.join("block.csv")).expect("block.csv");
    let labels: Vec<String> = reader
        .records()
        .map(|r| r.expect("block.csv record")[1].to_owned())
        .collect();
    let mut levels = labels.clone();
    levels.sort();
    levels.dedup();
    labels
        .iter()
        .map(|l| levels.binary_search(l).unwrap())
        .collect()
}

fn check(report: &mut Report, tol: Tol, what: &str, got: &[f64], want: &Matrix) {
    assert_eq!(got.len(), want.data.len(), "{what}: dimension mismatch");
    for j in 0..want.ncol {
        for i in 0..want.nrow {
            let (g, w) = (got[j * want.nrow + i], want.get(i, j));
            if g.is_nan() && w.is_nan() {
                continue;
            }
            report.check(
                tol,
                format_args!("{what}[{},{}]", want.row_names[i], want.col_names[j]),
                g,
                w,
            );
        }
    }
}

#[test]
fn duplicate_correlation_matches_golden() {
    for case in CASES {
        let dir = case_dir(case);
        let exprs = Matrix::read(&dir, "input_log2.csv");
        let design = Matrix::read(&dir, "design.csv");
        let block = read_block(&dir);
        let dc = duplicate_correlation(
            &exprs.data,
            exprs.nrow,
            exprs.ncol,
            &design.data,
            design.ncol,
            &block,
            0.15,
        )
        .expect("duplicateCorrelation");
        assert!(dc.warning.is_none(), "{case}: {:?}", dc.warning);
        let mut report = Report::new(&format!("{case}/dupcor"));
        let atanh = Matrix::read(&dir, "dupcor.csv");
        check(
            &mut report,
            REL_1E8,
            "atanh.correlations",
            &dc.atanh_correlations,
            &atanh,
        );
        let consensus = Matrix::read(&dir, "dupcor_scalars.csv");
        check(
            &mut report,
            REL_1E8,
            "consensus.correlation",
            &[dc.consensus_correlation],
            &consensus,
        );
        report.finish();
    }
}

#[test]
fn gls_fit_matches_golden() {
    for case in CASES {
        let dir = case_dir(case);
        let exprs = Matrix::read(&dir, "input_log2.csv");
        let design = Matrix::read(&dir, "design.csv");
        let block = read_block(&dir);
        let correlation = Matrix::read(&dir, "dupcor_scalars.csv").get(0, 0);
        let fit = lm_fit_matrix_block(
            &exprs.data,
            exprs.nrow,
            exprs.ncol,
            &design.data,
            design.ncol,
            &block,
            correlation,
        )
        .expect("lmFit(block)");
        let mut report = Report::new(&format!("{case}/lmfit"));
        let coef = Matrix::read(&dir, "lmfit_coefficients.csv");
        check(
            &mut report,
            ABS_1E10,
            "coefficients",
            &fit.coefficients,
            &coef,
        );
        let sd = Matrix::read(&dir, "lmfit_stdev_unscaled.csv");
        check(
            &mut report,
            REL_1E8,
            "stdev.unscaled",
            &fit.stdev_unscaled,
            &sd,
        );
        let scalars = Matrix::read(&dir, "lmfit_scalars.csv");
        for (name, got) in [
            ("sigma", &fit.sigma),
            ("df_residual", &fit.df_residual),
            ("Amean", &fit.amean),
        ] {
            for (i, (&g, &w)) in got.iter().zip(scalars.col(name)).enumerate() {
                if !(g.is_nan() && w.is_nan()) {
                    report.check(REL_1E8, format_args!("{name}[{}]", i + 1), g, w);
                }
            }
        }
        let cov = Matrix::read(&dir, "lmfit_cov_coefficients.csv");
        assert_eq!(cov.nrow, fit.rank, "{case}: rank");
        check(
            &mut report,
            REL_1E8,
            "cov.coefficients",
            &fit.cov_coefficients,
            &cov,
        );
        report.finish();
    }
}
