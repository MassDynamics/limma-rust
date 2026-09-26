//! Line-for-line port of R's `src/library/stats/src/lowess.c` (`clowess`, `lowest`) with
//! the `lowess()` R-level wrapper semantics (sort by x, `delta = 0.01 * diff(range(x))`,
//! `iter = 3`). See `tests/scalar_lowess.rs`.
