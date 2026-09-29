//! Pins the one deliberate deviation from limma 3.66.0 (README, "What the claim covers").
//!
//! `fitFDistUnequalDF1` with exactly two informative variances sets `prior.weights <- NULL` after
//! `PriorWeights` was already computed as TRUE, so R evaluates `w <- w * NULL`, which is
//! `numeric(0)`, and returns NaN scale and df2. The port drops the weights instead and fits the
//! two values unweighted. Reached only through the robust second pass (`ebayes.rs`, the call
//! with `Some(&fdr)`); the 2026-09-28 stats review judged it unreachable at production
//! feature counts.

use limma_core::ebayes::fit_f_dist_unequal_df1;

#[test]
fn two_informative_with_prior_weights_drops_the_weights() {
    let x = [0.8, 1.7, 0.0, 0.0];
    let df1 = [3.0, 5.0, 4.0, 4.0];
    let weighted =
        fit_f_dist_unequal_df1(&x, &df1, None, None, true, Some(&[0.2, 0.9, 1.0, 1.0])).unwrap();
    let plain = fit_f_dist_unequal_df1(&x, &df1, None, None, false, None).unwrap();

    // R: scale and df2 are NaN here. The port returns the unweighted, non-robust fit.
    assert_eq!(weighted.scale.len(), 1);
    assert!(weighted.scale[0].is_finite() && weighted.scale[0] > 0.0);
    assert!(!weighted.df2.is_nan());
    assert_eq!(weighted.scale, plain.scale);
    assert_eq!(weighted.df2.to_bits(), plain.df2.to_bits());
    assert!(weighted.df2_outlier.is_none() && weighted.df2_shrunk.is_none());
}

#[test]
fn fewer_than_two_informative_is_nan_as_in_r() {
    let fit = fit_f_dist_unequal_df1(&[0.8, 0.0, 0.0], &[3.0], None, None, false, None).unwrap();
    assert!(fit.scale[0].is_nan() && fit.df2.is_nan());
}
