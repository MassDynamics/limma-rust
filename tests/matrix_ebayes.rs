//! `eBayes` + `topTable` + `decideTests` against the `matrix/<case>/ebayes_rob*_trend*_*`
//! goldens, for all four robust/trend combinations, at the manifest tolerances: t/F/s2 rel
//! 1e-8, p rel 1e-8 with a 1e-12 floor, lods rel 1e-8, df/s2 prior rel 1e-6, logFC and CI abs
//! 1e-10, decideTests exact. The topTable goldens were written with `sort.by = "none"`, so
//! rows are matched by gene name rather than by position. The block goldens' downstream files
//! (`tests/golden/block/<case>/`) are held to the same bar.

mod common;

use common::matrix::{matrix_dir, Matrix, MATRIX_CASES};
use common::{Report, Tol};
use limma_core::block::lm_fit_matrix_block;
use limma_core::contrasts::contrasts_fit;
use limma_core::ebayes::{ebayes, EBayes, EBayesOptions};
use limma_core::fit::{lm_fit_matrix, MArrayLm};
use limma_core::toptable::{decide_tests_separate, top_table_f, top_table_t};

const ABS_1E10: Tol = Tol {
    rel: 0.0,
    abs_floor: 1e-10,
};
const REL_1E8: Tol = Tol {
    rel: 1e-8,
    abs_floor: 0.0,
};
const P_TOL: Tol = Tol {
    rel: 1e-8,
    abs_floor: 1e-12,
};
/// lods / B on the block cases: rel 1e-8 plus the logFC floor. lods is `log(p/(1-p))` (about
/// -4.6) plus terms of order one, so a value near zero carries ~1e-11 absolute rounding from
/// its inputs and fails a pure relative bar (bojkova_pairs gene 1: -7.6e-4, off by 1.1e-11).
const LODS_BLOCK: Tol = Tol {
    rel: 1e-8,
    abs_floor: 1e-10,
};
const REL_1E6: Tol = Tol {
    rel: 1e-6,
    abs_floor: 0.0,
};
const EXACT: Tol = Tol {
    rel: 0.0,
    abs_floor: 0.0,
};

fn load_case(dir: &std::path::Path) -> (Matrix, MArrayLm) {
    let exprs = Matrix::read(dir, "input_log2.csv");
    let design = Matrix::read(dir, "design.csv");
    let contrasts = Matrix::read(dir, "contrasts.csv");
    let fit = lm_fit_matrix(
        &exprs.data,
        exprs.nrow,
        exprs.ncol,
        &design.data,
        design.ncol,
    )
    .expect("lmFit");
    let cfit = contrasts_fit(&fit, &contrasts.data, contrasts.ncol).expect("contrasts.fit");
    (exprs, cfit)
}

fn check_matrix(report: &mut Report, tol: Tol, what: &str, got: &[f64], want: &Matrix) {
    assert_eq!(got.len(), want.data.len(), "{what}: dimension mismatch");
    for j in 0..want.ncol {
        for i in 0..want.nrow {
            let g = got[j * want.nrow + i];
            let w = want.get(i, j);
            if g.is_nan() && w.is_nan() {
                continue;
            }
            report.check(
                tol,
                format!("{what}[{},{}]", want.row_names[i], want.col_names[j]),
                g,
                w,
            );
        }
    }
}

fn check_vec(report: &mut Report, tol: Tol, what: &str, got: &[f64], want: &[f64]) {
    assert_eq!(got.len(), want.len(), "{what}: length mismatch");
    for (i, (&g, &w)) in got.iter().zip(want).enumerate() {
        if g.is_nan() && w.is_nan() {
            continue;
        }
        report.check(tol, format!("{what}[{}]", i + 1), g, w);
    }
}

fn index_by_name(m: &Matrix) -> std::collections::HashMap<&str, usize> {
    m.row_names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect()
}

fn recycle(v: &[f64], n: usize) -> Vec<f64> {
    if v.len() == 1 {
        vec![v[0]; n]
    } else {
        v.to_vec()
    }
}

fn check_ebayes(
    report: &mut Report,
    dir: &std::path::Path,
    prefix: &str,
    cfit: &MArrayLm,
    eb: &EBayes,
    lods_tol: Tol,
) {
    let n = cfit.ngenes;
    let t = Matrix::read(dir, &format!("{prefix}_t.csv"));
    check_matrix(report, REL_1E8, "t", &eb.t, &t);
    let p = Matrix::read(dir, &format!("{prefix}_p.csv"));
    check_matrix(report, P_TOL, "p.value", &eb.p_value, &p);
    let lods = Matrix::read(dir, &format!("{prefix}_lods.csv"));
    check_matrix(report, lods_tol, "lods", &eb.lods, &lods);
    let s = Matrix::read(dir, &format!("{prefix}_scalars.csv"));
    check_vec(report, REL_1E8, "s2.post", &eb.s2_post, s.col("s2_post"));
    check_vec(report, REL_1E6, "df.total", &eb.df_total, s.col("df_total"));
    check_vec(report, REL_1E8, "F", eb.f.as_ref().expect("F"), s.col("F"));
    check_vec(
        report,
        P_TOL,
        "F.p.value",
        eb.f_p_value.as_ref().expect("F.p"),
        s.col("F_p"),
    );
    check_vec(
        report,
        REL_1E6,
        "s2.prior",
        &recycle(&eb.s2_prior, n),
        s.col("s2_prior"),
    );
    check_vec(
        report,
        REL_1E6,
        "df.prior",
        &recycle(&eb.df_prior, n),
        s.col("df_prior"),
    );
}

fn check_toptables(
    report: &mut Report,
    dir: &std::path::Path,
    prefix: &str,
    exprs: &Matrix,
    cfit: &MArrayLm,
    eb: &EBayes,
    lods_tol: Tol,
) {
    for k in 0..cfit.ncoef {
        let want = Matrix::read(dir, &format!("{prefix}_toptable_{}.csv", k + 1));
        let rows = top_table_t(cfit, eb, k, 0.95).expect("topTable");
        assert_eq!(rows.len(), want.nrow);
        let by_name = index_by_name(&want);
        for row in &rows {
            let name = &exprs.row_names[row.gene];
            let i = by_name[name.as_str()];
            let tag = format!("toptable_{}[{name}]", k + 1);
            report.check(ABS_1E10, format!("{tag}.logFC"), row.log_fc, want.get(i, 0));
            report.check(ABS_1E10, format!("{tag}.CI.L"), row.ci_l, want.get(i, 1));
            report.check(ABS_1E10, format!("{tag}.CI.R"), row.ci_r, want.get(i, 2));
            report.check(
                REL_1E8,
                format!("{tag}.AveExpr"),
                row.ave_expr,
                want.get(i, 3),
            );
            report.check(REL_1E8, format!("{tag}.t"), row.t, want.get(i, 4));
            report.check(P_TOL, format!("{tag}.P.Value"), row.p_value, want.get(i, 5));
            report.check(
                P_TOL,
                format!("{tag}.adj.P.Val"),
                row.adj_p_value,
                want.get(i, 6),
            );
            report.check(lods_tol, format!("{tag}.B"), row.b, want.get(i, 7));
        }
    }
    let want = Matrix::read(dir, &format!("{prefix}_toptable_F.csv"));
    let rows = top_table_f(cfit, eb).expect("topTable F");
    assert_eq!(rows.len(), want.nrow);
    let nc = cfit.ncoef;
    assert_eq!(want.ncol, nc + 4);
    let by_name = index_by_name(&want);
    for row in &rows {
        let name = &exprs.row_names[row.gene];
        let i = by_name[name.as_str()];
        let tag = format!("toptable_F[{name}]");
        for j in 0..nc {
            report.check(
                ABS_1E10,
                format!("{tag}.coef{}", j + 1),
                row.coefs[j],
                want.get(i, j),
            );
        }
        report.check(
            REL_1E8,
            format!("{tag}.AveExpr"),
            row.ave_expr,
            want.get(i, nc),
        );
        report.check(REL_1E8, format!("{tag}.F"), row.f, want.get(i, nc + 1));
        report.check(
            P_TOL,
            format!("{tag}.P.Value"),
            row.p_value,
            want.get(i, nc + 2),
        );
        report.check(
            P_TOL,
            format!("{tag}.adj.P.Val"),
            row.adj_p_value,
            want.get(i, nc + 3),
        );
    }
}

fn check_grid(case: &str, dir: &std::path::Path, exprs: &Matrix, cfit: &MArrayLm, lods_tol: Tol) {
    for robust in [false, true] {
        for trend in [false, true] {
            let prefix = format!(
                "ebayes_rob{}_trend{}",
                if robust { "TRUE" } else { "FALSE" },
                if trend { "TRUE" } else { "FALSE" }
            );
            let mut report = Report::new(&format!("{case}/{prefix}"));
            let opts = EBayesOptions {
                robust,
                trend,
                ..EBayesOptions::default()
            };
            let eb = ebayes(cfit, &opts).expect("eBayes");
            check_ebayes(&mut report, dir, &prefix, cfit, &eb, lods_tol);
            check_toptables(&mut report, dir, &prefix, exprs, cfit, &eb, lods_tol);
            let want = Matrix::read(dir, &format!("{prefix}_decidetests.csv"));
            let got = decide_tests_separate(cfit, &eb, 0.05);
            check_matrix(&mut report, EXACT, "decideTests", &got, &want);
            report.finish();
        }
    }
}

#[test]
fn ebayes_matches_golden() {
    for case in MATRIX_CASES {
        let Some(dir) = matrix_dir("ebayes_matches_golden", case) else {
            return;
        };
        let (exprs, cfit) = load_case(&dir);
        check_grid(case, &dir, &exprs, &cfit, REL_1E8);
    }
}

/// The block cases in `tests/golden/block/` (see `block_goldens.R` there): genes with no
/// observation dropped, `lmFit(block = , correlation = )` at R's consensus correlation, then
/// the same contrasts / eBayes / topTable checks.
fn load_block_case(dir: &std::path::Path) -> (Matrix, MArrayLm) {
    let all = Matrix::read(dir, "input_log2.csv");
    let keep: Vec<usize> = (0..all.nrow)
        .filter(|&i| (0..all.ncol).any(|j| !all.get(i, j).is_nan()))
        .collect();
    let mut data = Vec::with_capacity(keep.len() * all.ncol);
    for j in 0..all.ncol {
        data.extend(keep.iter().map(|&i| all.get(i, j)));
    }
    let exprs = Matrix {
        file: all.file.clone(),
        row_names: keep.iter().map(|&i| all.row_names[i].clone()).collect(),
        col_names: all.col_names.clone(),
        nrow: keep.len(),
        ncol: all.ncol,
        data,
    };
    let design = Matrix::read(dir, "design.csv");
    let contrasts = Matrix::read(dir, "contrasts.csv");
    let correlation = Matrix::read(dir, "dupcor_scalars.csv").get(0, 0);
    let mut reader = csv::Reader::from_path(dir.join("block.csv")).expect("block.csv");
    let labels: Vec<String> = reader
        .records()
        .map(|r| r.expect("block.csv record")[1].to_owned())
        .collect();
    let mut levels = labels.clone();
    levels.sort();
    levels.dedup();
    let block: Vec<usize> = labels
        .iter()
        .map(|l| levels.binary_search(l).unwrap())
        .collect();
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
    let cfit = contrasts_fit(&fit, &contrasts.data, contrasts.ncol).expect("contrasts.fit");
    (exprs, cfit)
}

#[test]
fn block_ebayes_matches_golden() {
    for case in [
        "techrep_synth",
        "techrep_synth_na",
        "bojkova_pairs",
        "techrep_mixed",
    ] {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/golden/block")
            .join(case);
        let (exprs, cfit) = load_block_case(&dir);
        let mut report = Report::new(&format!("{case}/contrasts"));
        for (what, got, file) in [
            (
                "coefficients",
                &cfit.coefficients,
                "contrasts_coefficients.csv",
            ),
            (
                "stdev.unscaled",
                &cfit.stdev_unscaled,
                "contrasts_stdev_unscaled.csv",
            ),
        ] {
            let tol = if what == "coefficients" {
                ABS_1E10
            } else {
                REL_1E8
            };
            check_matrix(&mut report, tol, what, got, &Matrix::read(&dir, file));
        }
        report.finish();
        check_grid(case, &dir, &exprs, &cfit, LODS_BLOCK);
    }
}
