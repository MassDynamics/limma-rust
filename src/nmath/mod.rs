//! R's distribution functions and gamma family, ported from `R-4.5.0/src/nmath`.
//!
//! Every function here is a port of a named R C function, not an equivalent: R's tails are
//! the spec, and the scalar golden corpus is what proves it (see `tests/scalar_goldens.rs`).
//! Each public function's doc comment names the R source file and function it comes from.
//!
//! The two modules below carry R's private headers, and stay `pub(crate)` for the same
//! reason those headers are private:
//!
//! - `arith`  `fmax2`, `fmin2`, `R_pow`, `R_pow_di` (`nmath.h` / `arithmetic.c`).
//! - `consts` the `M_*` constants `nmath.h` pulls in from `Rmath.h`.
//! - `dpq`    the `p`/`q` helpers of `dpq.h`, plus `nmath.h`'s `R_forceint`.

// R's private headers land a card ahead of the functions that call them, so the whole of
// `dpq` and most of `consts` has no caller yet. Delete this once the ports are in — after
// that, an unused helper here is a real finding.
#![allow(dead_code, unused_imports, unused_macros)]

pub(crate) mod arith;
pub(crate) mod consts;
pub(crate) mod dpq;

pub mod f;
pub mod gamma;
pub mod pbeta;
pub mod pgamma;
pub mod pnorm;
pub mod qbeta;
pub mod t;
pub mod toms708;
pub use pnorm::{dnorm, pnorm, qnorm};
