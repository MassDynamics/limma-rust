//! The block path of `lmFit`: `duplicateCorrelation(object, design, block = )` (limma 3.68.5,
//! `R/dups.R`) with statmod 1.5.2's `mixedModel2Fit(only.varcomp = TRUE)` and `glmgam.fit`,
//! and `gls.series(M, design, block = , correlation = )` (`R/lmfit.R`). No weights and no
//! duplicate spots (`ndups`): only the block arm of each function.
//!
//! `block` is a vector of integer codes, one per array. R's `factor(block)` orders the levels,
//! so callers must code blocks in the order R would sort their labels (the only place the order
//! shows is which level `duplicateCorrelation` drops when it checks whether the design already
//! encodes the block).
//!
//! `mixedModel2Fit` takes `La.svd(QtZ, nu = mq)`; here the squared singular values and left
//! singular vectors come from the symmetric eigen-decomposition of `QtZ %*% t(QtZ)`. The vectors
//! of a repeated value are only defined up to a rotation, but the gamma fit only sees the sum of
//! `dy` over a group of equal rows of `dx`, which is rotation-invariant, so the two agree to
//! rounding rather than bit for bit.

use crate::fit::{row_means_na_rm, MArrayLm, LM_TOL};
use crate::linalg::{chol_upper, eigen_symmetric, matmul, mean, mean_trim};
use crate::linpack::{lm_fit, qr_decompose};
use crate::{LimmaError, Result};

/// `duplicateCorrelation(...)$consensus.correlation` and `$atanh.correlations`.
#[derive(Debug, Clone)]
pub struct DupCor {
    pub consensus_correlation: f64,
    /// One per gene; `NaN` where the per-gene fit was skipped or failed.
    pub atanh_correlations: Vec<f64>,
    /// The warning R raises when it returns a zero correlation without fitting.
    pub warning: Option<&'static str>,
}

/// `factor(block)`'s sorted levels and each array's 0-based level index.
fn block_levels(block: &[usize]) -> (Vec<usize>, Vec<usize>) {
    let mut levels = block.to_vec();
    levels.sort_unstable();
    levels.dedup();
    let idx = block
        .iter()
        .map(|b| levels.binary_search(b).expect("level"))
        .collect();
    (levels, idx)
}

/// `duplicateCorrelation(M, design, block = block, trim = trim)`. `m` is `ngenes x narrays` and
/// `design` is `narrays x nbeta`, both column-major; `NaN` plays `NA`.
pub fn duplicate_correlation(
    m: &[f64],
    ngenes: usize,
    narrays: usize,
    design: &[f64],
    nbeta: usize,
    block: &[usize],
    trim: f64,
) -> Result<DupCor> {
    if m.len() != ngenes * narrays {
        return Err(LimmaError::Invalid("exprs has the wrong length".into()));
    }
    if design.len() != narrays * nbeta {
        return Err(LimmaError::Invalid(
            "Number of rows of design matrix does not match number of arrays".into(),
        ));
    }
    if block.len() != narrays {
        return Err(LimmaError::Invalid(
            "Length of block does not match number of arrays".into(),
        ));
    }
    let zero = |warning| DupCor {
        consensus_correlation: 0.0,
        atanh_correlations: vec![0.0; ngenes],
        warning: Some(warning),
    };

    let (levels, level_of) = block_levels(block);
    let mut sizes = vec![0usize; levels.len()];
    for &l in &level_of {
        sizes[l] += 1;
    }
    let max_block_size = sizes.iter().copied().max().unwrap_or(0);
    if max_block_size == 1 {
        return Ok(zero(
            "Blocks all of size 1: setting intrablock correlation to zero.",
        ));
    }
    // Block factor already in the design? model.matrix(~factor(block))[, -1] against qr(design).
    let qr = qr_decompose(design, narrays, nbeta, LM_TOL);
    let mut max_abs = f64::NEG_INFINITY;
    for l in 1..levels.len() {
        let col: Vec<f64> = level_of
            .iter()
            .map(|&k| if k == l { 1.0 } else { 0.0 })
            .collect();
        for v in &qr.qty(&col)[qr.rank..] {
            max_abs = max_abs.max(v.abs());
        }
    }
    if max_abs < 1e-8 {
        return Ok(zero(
            "Block factor already encoded in the design matrix: setting intrablock correlation to zero.",
        ));
    }

    let mut rho = vec![f64::NAN; ngenes];
    let mut x = Vec::new();
    let mut z = Vec::new();
    for (g, rho_g) in rho.iter_mut().enumerate() {
        let obs: Vec<usize> = (0..narrays)
            .filter(|&j| m[j * ngenes + g].is_finite())
            .collect();
        let nobs = obs.len();
        // Levels of factor(block[o]), in sorted order.
        let mut present: Vec<usize> = obs.iter().map(|&j| level_of[j]).collect();
        present.sort_unstable();
        present.dedup();
        let nblocks = present.len();
        if !(nobs > nbeta + 2 && nblocks > 1 && nblocks + 1 < nobs) {
            continue;
        }
        let y: Vec<f64> = obs.iter().map(|&j| m[j * ngenes + g]).collect();
        x.clear();
        for k in 0..nbeta {
            x.extend(obs.iter().map(|&j| design[k * narrays + j]));
        }
        z.clear();
        for &l in &present {
            z.extend(
                obs.iter()
                    .map(|&j| if level_of[j] == l { 1.0 } else { 0.0 }),
            );
        }
        // tryCatch(..., error = function(e) NA)
        if let Ok(s) = mixed_model2_fit_varcomp(&y, &x, nbeta, &z, nblocks, 20) {
            if !s[0].is_nan() {
                *rho_g = s[1] / (s[0] + s[1]);
            }
        }
    }

    let rhomax = 0.99;
    let rhomin = 1.0 / (1.0 - max_block_size as f64) + 0.01;
    let finite = || rho.iter().copied().filter(|r| !r.is_nan());
    if finite().fold(0.0, f64::min) < rhomin {
        for r in rho.iter_mut().filter(|r| **r < rhomin) {
            *r = rhomin;
        }
    }
    let finite = || rho.iter().copied().filter(|r| !r.is_nan());
    if finite().fold(0.0, f64::max) > rhomax {
        for r in rho.iter_mut().filter(|r| **r > rhomax) {
            *r = rhomax;
        }
    }
    let arho: Vec<f64> = rho.iter().map(|r| r.atanh()).collect();
    let kept: Vec<f64> = arho.iter().copied().filter(|a| !a.is_nan()).collect();
    Ok(DupCor {
        consensus_correlation: mean_trim(&kept, trim).tanh(),
        atanh_correlations: arho,
        warning: None,
    })
}

/// `statmod::mixedModel2Fit(y, X, Z, only.varcomp = TRUE, maxit = maxit)$varcomp`:
/// `(Residual, Block)`. `x` is `n x nx`, `z` is `n x nz`, both column-major. An `Err` is what R
/// would raise (and `duplicateCorrelation` turns into `NA`).
pub fn mixed_model2_fit_varcomp(
    y: &[f64],
    x: &[f64],
    nx: usize,
    z: &[f64],
    nz: usize,
    maxit: usize,
) -> Result<[f64; 2]> {
    let mx = y.len();
    let qr = qr_decompose(x, mx, nx, LM_TOL);
    let r = qr.rank;
    let mq = mx - r;
    if mq == 0 {
        return Ok([f64::NAN, f64::NAN]);
    }
    // QtZ = effects[(r+1):mx, 1:nz] of lm.fit(X, cbind(Z, y)), and the y column's tail.
    let mut qtz = vec![0.0; mq * nz];
    for k in 0..nz {
        let e = qr.qty(&z[k * mx..(k + 1) * mx]);
        qtz[k * mq..(k + 1) * mq].copy_from_slice(&e[r..]);
    }
    let qy = qr.qty(y)[r..].to_vec();

    let (d, dy) = svd_d2_dy(&qtz, mq, nz, &qy);
    let mut dx = vec![1.0; mq];
    dx.extend_from_slice(&d);

    if dy.iter().any(|v| !v.is_finite()) {
        return Err(LimmaError::Invalid("NA/NaN/Inf in 'y'".into()));
    }
    let dfit = lm_fit(&dx, mq, 2, &dy, LM_TOL)?;
    let mut varcomp = [dfit.coefficients[0], dfit.coefficients[1]];
    let nonzero = d.iter().filter(|v| v.abs() > 1e-15).count();
    if mq > 2 && nonzero > 1 && var(&d) > 1e-15 {
        let start = if dfit.fitted_values.iter().all(|&f| f >= 0.0) {
            varcomp.to_vec()
        } else {
            vec![mean(&dy), 0.0]
        };
        let beta = glmgam_fit(&dx, mq, 2, &dy, &start, 1e-6, maxit)?;
        varcomp = [beta[0], beta[1]];
    }
    Ok(varcomp)
}

/// `s <- La.svd(QtZ, nu = mq); d <- s$d^2; dy <- drop(crossprod(s$u, Qy))^2`, padded with zeros
/// to `mq` as `mixedModel2Fit` does. Squared singular values of R's structural zeros are near
/// 1e-32; eigenvalues carry noise of order `mq * eps * max`, which could pass the later
/// `abs(d) > 1e-15` gate, so anything at noise level is set to zero.
///
/// With fewer blocks than residual dimensions the work is the `nz x nz` problem `t(QtZ) %*% QtZ`:
/// `u = QtZ v / sqrt(lambda)` for the nonzero directions, and the rest of `Qy` lies in the null
/// space, where every row of `dx` is `(1, 0)`. The fits downstream depend on `dy` only through its
/// sum over equal rows of `dx` (the null basis is arbitrary in R too), so that mass is spread
/// evenly over the null rows.
fn svd_d2_dy(qtz: &[f64], mq: usize, nz: usize, qy: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let noise = |top: f64| mq as f64 * f64::EPSILON * top.abs();
    if nz >= mq {
        let mut zzt = vec![0.0; mq * mq];
        for j in 0..mq {
            for i in 0..mq {
                let mut s = 0.0;
                for k in 0..nz {
                    s += qtz[k * mq + i] * qtz[k * mq + j];
                }
                zzt[j * mq + i] = s;
            }
        }
        let (values, u) = eigen_symmetric(&zzt, mq);
        let floor = noise(values[0]);
        let d = values
            .iter()
            .map(|&v| if v > floor { v } else { 0.0 })
            .collect();
        let dy = (0..mq)
            .map(|j| {
                let s: f64 = (0..mq).map(|i| u[j * mq + i] * qy[i]).sum();
                s * s
            })
            .collect();
        return (d, dy);
    }
    let mut ztz = vec![0.0; nz * nz];
    for b in 0..nz {
        for a in 0..nz {
            ztz[b * nz + a] = (0..mq).map(|i| qtz[a * mq + i] * qtz[b * mq + i]).sum();
        }
    }
    let (values, v) = eigen_symmetric(&ztz, nz);
    let floor = noise(values[0]);
    let mut d = vec![0.0; mq];
    let mut dy = vec![0.0; mq];
    let mut resid = qy.to_vec();
    let mut k = 0;
    for (j, &lambda) in values.iter().enumerate() {
        if lambda <= floor {
            continue;
        }
        let sigma = lambda.sqrt();
        let u: Vec<f64> = (0..mq)
            .map(|i| {
                (0..nz)
                    .map(|b| qtz[b * mq + i] * v[j * nz + b])
                    .sum::<f64>()
                    / sigma
            })
            .collect();
        let c: f64 = u.iter().zip(qy).map(|(a, b)| a * b).sum();
        for (r, a) in resid.iter_mut().zip(&u) {
            *r -= c * a;
        }
        d[k] = lambda;
        dy[k] = c * c;
        k += 1;
    }
    let null = resid.iter().map(|r| r * r).sum::<f64>() / (mq - k) as f64;
    for v in &mut dy[k..] {
        *v = null;
    }
    (d, dy)
}

/// R's `var(x)` for a plain vector.
fn var(x: &[f64]) -> f64 {
    let n = x.len();
    if n < 2 {
        return f64::NAN;
    }
    let m = mean(x);
    x.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / (n - 1) as f64
}

/// R's `if (cond)`: an `NA` condition is an error.
fn cond(c: Option<bool>) -> Result<bool> {
    c.ok_or_else(|| LimmaError::Invalid("missing value where TRUE/FALSE needed".into()))
}

/// `a <= b`, `a < b` with R's `NA` for a `NaN` operand.
fn le(a: f64, b: f64) -> Option<bool> {
    (!a.is_nan() && !b.is_nan()).then_some(a <= b)
}
fn lt(a: f64, b: f64) -> Option<bool> {
    (!a.is_nan() && !b.is_nan()).then_some(a < b)
}
/// R's `||`: `TRUE || NA` is `TRUE`, `FALSE || NA` is `NA`.
fn or(a: Option<bool>, b: impl FnOnce() -> Option<bool>) -> Option<bool> {
    match a {
        Some(true) => Some(true),
        Some(false) => b(),
        None => match b() {
            Some(true) => Some(true),
            _ => None,
        },
    }
}

/// `max(x)` with R's `NA` propagation.
fn rmax(x: &[f64]) -> f64 {
    let mut m = f64::NEG_INFINITY;
    for &v in x {
        if v.is_nan() {
            return f64::NAN;
        }
        m = m.max(v);
    }
    m
}

/// `deviance.gamma(y, mu)` inside `glmgam.fit`.
fn deviance_gamma(y: &[f64], mu: &[f64]) -> Result<f64> {
    if mu.iter().any(|v| v.is_nan()) && !mu.iter().any(|&v| v < 0.0) {
        return cond(None).map(|_| 0.0);
    }
    if mu.iter().any(|&v| v < 0.0) {
        return Ok(f64::INFINITY);
    }
    let mut s = 0.0;
    let mut any_kept = false;
    for (&yi, &mi) in y.iter().zip(mu) {
        if yi < 1e-15 && mi < 1e-15 {
            continue;
        }
        any_kept = true;
        s += (yi - mi) / mi - (yi / mi).ln();
    }
    Ok(if any_kept { 2.0 * s } else { 0.0 })
}

/// Solve `(R'R) x = b` for an upper-triangular `R` (`backsolve(R, backsolve(R, b, transpose =
/// TRUE))`).
fn chol_solve(r: &[f64], n: usize, b: &[f64]) -> Vec<f64> {
    let mut z = forward_solve_t(r, n, b);
    for i in (0..n).rev() {
        let mut t = z[i];
        for k in (i + 1)..n {
            t -= r[k * n + i] * z[k];
        }
        z[i] = t / r[i * n + i];
    }
    z
}

/// `backsolve(U, b, transpose = TRUE)`: solve `t(U) x = b` for upper-triangular `U` (`n x n`),
/// in BLAS `dtrsm`'s order.
fn forward_solve_t(u: &[f64], n: usize, b: &[f64]) -> Vec<f64> {
    let mut x = b.to_vec();
    for i in 0..n {
        let mut t = x[i];
        for k in 0..i {
            t -= u[i * n + k] * x[k];
        }
        x[i] = t / u[i * n + i];
    }
    x
}

/// `statmod::glmgam.fit(X, y, coef.start = start, tol, maxit)$coefficients`: a gamma GLM with
/// identity link by Levenberg-damped scoring. `x` is `n x p`, column-major.
pub fn glmgam_fit(
    x: &[f64],
    n: usize,
    p: usize,
    y: &[f64],
    start: &[f64],
    tol: f64,
    maxit: usize,
) -> Result<Vec<f64>> {
    if p > n {
        return Err(LimmaError::Invalid("More columns than rows in X".into()));
    }
    if n == 0 {
        return Ok(Vec::new());
    }
    if !(y.iter().all(|v| v.is_finite()) || x.iter().all(|v| v.is_finite())) {
        return Err(LimmaError::Invalid(
            "All values must be finite and non-missing".into(),
        ));
    }
    if y.iter().any(|&v| v < 0.0) {
        return Err(LimmaError::Invalid("y must be non-negative".into()));
    }
    if rmax(y) == 0.0 {
        return Ok(vec![0.0; p]);
    }
    let mut beta = start.to_vec();
    let mut mu = matmul(x, n, p, &beta, 1);
    if mu.iter().any(|&v| v < 0.0) {
        return Err(LimmaError::Invalid(
            "Starting values give negative fitted values".into(),
        ));
    }
    let mut dev = deviance_gamma(y, &mu)?;
    let mut iter = 0;
    let mut lambda = 0.0;
    loop {
        iter += 1;
        let mut v: Vec<f64> = mu.iter().map(|m| m * m).collect();
        let vfloor = rmax(&v) / 1e3;
        for vi in v.iter_mut() {
            // pmax: an NA stays NA
            if !vi.is_nan() && *vi < vfloor {
                *vi = vfloor;
            }
        }
        let mut xvx = vec![0.0; p * p];
        for j in 0..p {
            for i in 0..p {
                let mut s = 0.0;
                for k in 0..n {
                    s += x[i * n + k] * (x[j * n + k] / v[k]);
                }
                xvx[j * p + i] = s;
            }
        }
        let diag: Vec<f64> = (0..p).map(|i| xvx[i * p + i]).collect();
        let maxinfo = rmax(&diag);
        if iter == 1 {
            lambda = mean(&diag).abs() / p as f64;
        }
        let resid: Vec<f64> = (0..n).map(|k| (y[k] - mu[k]) / v[k]).collect();
        let dl: Vec<f64> = (0..p)
            .map(|i| (0..n).map(|k| x[i * n + k] * resid[k]).sum())
            .collect();
        let betaold = beta.clone();
        let devold = dev;
        let mut lev = 0;
        let mut dbeta;
        loop {
            lev += 1;
            let mut a = xvx.clone();
            for i in 0..p {
                a[i * p + i] += lambda;
            }
            let r = chol_upper(&a, p)?;
            dbeta = chol_solve(&r, p, &dl);
            beta = betaold.iter().zip(&dbeta).map(|(b, d)| b + d).collect();
            mu = matmul(x, n, p, &beta, 1);
            dev = deviance_gamma(y, &mu)?;
            if cond(or(le(dev, devold), || lt(dev / rmax(&mu), 1e-15)))? {
                break;
            }
            if cond(lt(1e15, lambda / maxinfo))? {
                beta = betaold.clone();
                break;
            }
            lambda *= 2.0;
        }
        if cond(lt(1e15, lambda / maxinfo))? {
            break;
        }
        if lev == 1 {
            lambda /= 10.0;
        }
        let step: f64 = dl.iter().zip(&dbeta).map(|(a, b)| a * b).sum();
        if cond(or(lt(step, tol), || lt(dev / rmax(&mu), 1e-15)))? {
            break;
        }
        if iter > maxit {
            break;
        }
    }
    Ok(beta)
}

/// `gls.series(M, design, block = block, correlation = correlation)`: generalised least squares
/// with a common intra-block correlation. One whitened `lm.fit` for all genes when every value
/// is finite, otherwise one fit per gene on its finite observations.
pub fn gls_series(
    m: &[f64],
    ngenes: usize,
    narrays: usize,
    design: &[f64],
    nbeta: usize,
    block: &[usize],
    correlation: f64,
) -> Result<MArrayLm> {
    let n = narrays;
    let p = nbeta;
    if design.len() != n * p {
        return Err(LimmaError::Invalid(
            "Number of rows of design matrix does not match number of arrays".into(),
        ));
    }
    if block.len() != n {
        return Err(LimmaError::Invalid(
            "Length of block does not match number of arrays".into(),
        ));
    }
    if correlation.is_nan() {
        return Err(LimmaError::Invalid(
            "missing value where TRUE/FALSE needed".into(),
        ));
    }
    if correlation.abs() >= 1.0 {
        return Err(LimmaError::Invalid(
            "correlation is 1 or -1, so the model is degenerate".into(),
        ));
    }
    let mut cormatrix = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            cormatrix[j * n + i] = if i == j {
                1.0
            } else if block[i] == block[j] {
                correlation
            } else {
                0.0
            };
        }
    }

    let mut coefficients = vec![f64::NAN; ngenes * p];
    let mut stdev_unscaled = vec![f64::NAN; ngenes * p];
    let mut sigma = vec![f64::NAN; ngenes];
    let mut df_residual = vec![0.0; ngenes];

    let chol_v = chol_upper(&cormatrix, n)?;
    let whiten = |u: &[f64], k: usize, cols: &[f64], ncols: usize| -> Vec<f64> {
        let mut out = Vec::with_capacity(k * ncols);
        for c in 0..ncols {
            out.extend(forward_solve_t(u, k, &cols[c * k..(c + 1) * k]));
        }
        out
    };
    let xw = whiten(&chol_v, n, design, p);
    let qr = qr_decompose(&xw, n, p, LM_TOL);
    let rank = qr.rank;
    let cov_coefficients = qr.chol2inv();

    if m.iter().all(|v| v.is_finite()) {
        let df = n - rank;
        let mut y = vec![0.0; n];
        for g in 0..ngenes {
            for (j, yj) in y.iter_mut().enumerate() {
                *yj = m[j * ngenes + g];
            }
            let yw = forward_solve_t(&chol_v, n, &y);
            let qty = qr.qty(&yw);
            let coef = qr.coef_pivoted(&qty)?;
            for (k, &b) in coef.iter().enumerate() {
                coefficients[qr.pivot[k] * ngenes + g] = b;
            }
            if df > 0 {
                // colMeans(effects[-(1:rank), ]^2)
                let ss: f64 = qty[rank..].iter().map(|e| e * e).sum();
                sigma[g] = (ss / df as f64).sqrt();
            }
            df_residual[g] = df as f64;
        }
        for k in 0..rank {
            let sd = cov_coefficients[k * rank + k].sqrt();
            for g in 0..ngenes {
                stdev_unscaled[qr.pivot[k] * ngenes + g] = sd;
            }
        }
    } else {
        let mut x = Vec::with_capacity(n * p);
        let mut v = Vec::with_capacity(n * n);
        for g in 0..ngenes {
            let obs: Vec<usize> = (0..n).filter(|&j| m[j * ngenes + g].is_finite()).collect();
            let nobs = obs.len();
            if nobs == 0 {
                continue;
            }
            v.clear();
            for &j in &obs {
                v.extend(obs.iter().map(|&i| cormatrix[j * n + i]));
            }
            let u = chol_upper(&v, nobs)?;
            let y: Vec<f64> = obs.iter().map(|&j| m[j * ngenes + g]).collect();
            let yw = forward_solve_t(&u, nobs, &y);
            x.clear();
            for k in 0..p {
                x.extend(obs.iter().map(|&j| design[k * n + j]));
            }
            if x.iter().all(|&v| v == 0.0) {
                df_residual[g] = nobs as f64;
                let w = 1.0 / nobs as f64;
                sigma[g] = yw.iter().map(|e| w * (e * e)).sum::<f64>().sqrt();
                continue;
            }
            let xg = whiten(&u, nobs, &x, p);
            let out = lm_fit(&xg, nobs, p, &yw, LM_TOL)?;
            for k in 0..p {
                coefficients[k * ngenes + g] = out.coefficients[k];
            }
            let r = out.rank;
            let cov = out.qr.chol2inv();
            for k in 0..r {
                stdev_unscaled[out.qr.pivot[k] * ngenes + g] = cov[k * r + k].sqrt();
            }
            df_residual[g] = out.df_residual as f64;
            if out.df_residual > 0 {
                let w = 1.0 / out.df_residual as f64;
                sigma[g] = out
                    .residuals
                    .iter()
                    .map(|e| w * (e * e))
                    .sum::<f64>()
                    .sqrt();
            }
        }
    }

    Ok(MArrayLm {
        ngenes,
        ncoef: p,
        narrays: n,
        coefficients,
        stdev_unscaled,
        sigma,
        df_residual,
        cov_coefficients,
        pivot: qr.pivot,
        rank,
        amean: Vec::new(),
        design: design.to_vec(),
    })
}

/// `lmFit(exprs, design, block = block, correlation = correlation)`.
pub fn lm_fit_matrix_block(
    exprs: &[f64],
    ngenes: usize,
    narrays: usize,
    design: &[f64],
    ncoef: usize,
    block: &[usize],
    correlation: f64,
) -> Result<MArrayLm> {
    if exprs.len() != ngenes * narrays {
        return Err(LimmaError::Invalid("exprs has the wrong length".into()));
    }
    if design.len() != narrays * ncoef {
        return Err(LimmaError::Invalid(
            "row dimension of design doesn't match column dimension of data object".into(),
        ));
    }
    if design.iter().any(|v| !v.is_finite()) {
        return Err(LimmaError::Invalid(
            "design must not contain NA or Inf".into(),
        ));
    }
    let amean = row_means_na_rm(exprs, ngenes, narrays);
    let mut fit = gls_series(exprs, ngenes, narrays, design, ncoef, block, correlation)?;
    fit.amean = amean;
    Ok(fit)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_correlation_gls_is_ols() {
        let m = [1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 3.0, 6.0, 9.0, 4.0, 8.0, 12.5];
        let design = [1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 2.0, 3.0];
        let block = [0, 0, 1, 1];
        let ols = crate::fit::lm_series(&m, 3, 4, &design, 2).unwrap();
        let gls = gls_series(&m, 3, 4, &design, 2, &block, 0.0).unwrap();
        for (a, b) in ols.coefficients.iter().zip(&gls.coefficients) {
            assert!((a - b).abs() < 1e-12);
        }
        for (a, b) in ols.sigma.iter().zip(&gls.sigma) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn degenerate_correlation_is_an_error() {
        let m = [1.0, 2.0, 3.0, 4.0];
        let design = [1.0; 4];
        assert!(gls_series(&m, 1, 4, &design, 1, &[0, 0, 1, 1], 1.0).is_err());
    }

    #[test]
    fn singleton_blocks_give_zero_with_a_warning() {
        let m = [1.0, 2.0, 3.0, 4.0];
        let design = [1.0; 4];
        let dc = duplicate_correlation(&m, 1, 4, &design, 1, &[0, 1, 2, 3], 0.15).unwrap();
        assert_eq!(dc.consensus_correlation, 0.0);
        assert!(dc.warning.unwrap().starts_with("Blocks all of size 1"));
    }

    #[test]
    fn small_side_svd_matches_the_full_eigen_problem() {
        // mq = 7, nz = 3, third column = first + second so one squared singular value is zero.
        let (mq, nz) = (7, 3);
        let mut qtz: Vec<f64> = (0..mq * 2)
            .map(|i| ((i * 37 + 11) % 17) as f64 / 7.0 - 1.0)
            .collect();
        let third: Vec<f64> = (0..mq).map(|i| qtz[i] + qtz[mq + i]).collect();
        qtz.extend(third);
        let qy: Vec<f64> = (0..mq).map(|i| ((i * 13 + 5) % 11) as f64 - 4.5).collect();
        let (d, dy) = svd_d2_dy(&qtz, mq, nz, &qy);

        let zzt: Vec<f64> = (0..mq * mq)
            .map(|ji| {
                let (j, i) = (ji / mq, ji % mq);
                (0..nz).map(|k| qtz[k * mq + i] * qtz[k * mq + j]).sum()
            })
            .collect();
        let (values, u) = eigen_symmetric(&zzt, mq);
        let full_dy: Vec<f64> = (0..mq)
            .map(|j| (0..mq).map(|i| u[j * mq + i] * qy[i]).sum::<f64>().powi(2))
            .collect();
        assert_eq!(d.iter().filter(|&&v| v > 0.0).count(), 2);
        for j in 0..2 {
            assert!((d[j] - values[j]).abs() < 1e-12 * values[0]);
            assert!((dy[j] - full_dy[j]).abs() < 1e-12 * full_dy[j].max(1.0));
        }
        let null: f64 = dy[2..].iter().sum();
        let full_null: f64 = full_dy[2..].iter().sum();
        assert!((null - full_null).abs() < 1e-12 * full_null);
    }

    #[test]
    fn block_in_design_gives_zero_with_a_warning() {
        let m = [1.0, 2.0, 3.0, 4.0];
        let design = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0];
        let dc = duplicate_correlation(&m, 1, 4, &design, 2, &[0, 0, 1, 1], 0.15).unwrap();
        assert_eq!(dc.consensus_correlation, 0.0);
        assert!(dc
            .warning
            .unwrap()
            .starts_with("Block factor already encoded"));
    }
}
