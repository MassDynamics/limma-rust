//! `limma-core`: a pure-Rust port of the limma statistical core used by the Mass Dynamics
//! pairwise-comparison job. No R anywhere.
//!
//! The golden corpus produced by the R implementation IS the spec. Every module here is
//! checked against it (see `tests/` and the repo `README.md`), under the tolerance policy
//! recorded in each corpus manifest.
//!
//! Module map (filled in by the ordered forge cards, see `docs/cards.md`):
//!
//! - `nmath`     R's distribution functions and gamma family (`pt`, `qt`, `pf`, `qf`,
//!               `pnorm`, `qnorm`, `pbeta`, `qbeta`, `pgamma`, `qgamma`, `pchisq`, `qchisq`,
//!               `digamma`, `trigamma`) plus limma's `trigammaInverse`. Ported from R's
//!               `src/nmath`, not from a generic stats crate: tails must match R within the
//!               corpus tolerances.
//! - `lowess`    line-for-line port of R's C `clowess` (used by `eBayes(trend=TRUE)`).
//! - `lm`        `lmFit` / `lm.series` (QR via faer; NA-aware per-gene refits).
//! - `contrasts` `contrasts.fit`, incl. the non-orthogonal cov/Cholesky route.
//! - `ebayes`    `eBayes`, `squeezeVar`, `fitFDist`, `fitFDistRobustly`.
//! - `toptable`  `topTable` (CIs, BH adjust, moderated F), `decideTests`.
//! - `camera`    `camera` parametric path (used by the enrichment job).

pub const CRATE_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Errors surfaced to the Python layer. Messages that reach users must match the R
/// package's `md_error` texts where a corpus manifest records an `expected_error`.
#[derive(Debug, thiserror::Error)]
pub enum LimmaError {
    #[error("{0}")]
    Invalid(String),
}

pub type Result<T> = std::result::Result<T, LimmaError>;
