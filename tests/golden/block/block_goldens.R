# block_goldens.R
#
# Goldens for the block path: duplicateCorrelation(block = ) and lmFit(block = , correlation = )
# (gls.series). Writes one directory per case next to this script. Doubles are printed with
# %.17g so they round-trip exactly.
#
# Run from the limma-rust root in the corpus R image (R 4.5.0, limma 3.68.5, statmod 1.5.2):
#
#   docker run --rm --platform linux/amd64 -e LANG=C.UTF-8 -e LC_COLLATE=C -v $PWD:/w -w /w \
#     md-flexi-r45-limma368 Rscript tests/golden/block/block_goldens.R
#
# Cases:
#   techrep_synth     400 genes, 4 conditions x 3 biological samples x 2 technical replicates,
#                     all finite, design ~0 + condition. Block = biological sample.
#   techrep_synth_na  the same data with ~12% values missing at random plus a few genes too
#                     sparse for duplicateCorrelation to fit, design ~condition (intercept).
#   bojkova_pairs     the matrix/conditions_only slice (600 genes x 24 samples, public), its
#                     design, and blocks of two consecutive samples within each condition.
#   techrep_mixed     300 genes, 4 conditions x 4 biological samples with 3, 2, 3, 2 technical
#                     replicates (40 arrays, MaxBlockSize 3), design ~0 + condition. Gene 1 has
#                     near-identical replicates (rho clamped to rhomax 0.99); gene 2 is cut to
#                     pairs by NA with opposite-sign pair deviations (rho clamped to rhomin -0.49).

suppressPackageStartupMessages(library(limma))

out_root <- "tests/golden/block"

fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))

write_mat <- function(m, dir, file) {
  m <- as.matrix(m)
  rn <- rownames(m)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(m)))
  out <- data.frame(row = rn, apply(m, 2, fmt), check.names = FALSE, stringsAsFactors = FALSE)
  if (is.null(colnames(m))) colnames(out)[-1] <- paste0("V", seq_len(ncol(m)))
  write.csv(out, file.path(dir, file), row.names = FALSE, quote = 1)
}

write_case <- function(id, exprs, design, block) {
  dir <- file.path(out_root, id)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  rownames(exprs) <- as.character(seq_len(nrow(exprs)))
  write_mat(exprs, dir, "input_log2.csv")
  write_mat(design, dir, "design.csv")
  write.csv(data.frame(row = as.character(seq_along(block)), block = block),
            file.path(dir, "block.csv"), row.names = FALSE, quote = 1)

  dc <- duplicateCorrelation(exprs, design, block = block)
  write_mat(cbind(atanh_correlation = dc$atanh.correlations), dir, "dupcor.csv")
  write_mat(cbind(consensus_correlation = dc$consensus.correlation), dir, "dupcor_scalars.csv")

  fit <- lmFit(exprs, design, block = block, correlation = dc$consensus.correlation)
  write_mat(fit$coefficients, dir, "lmfit_coefficients.csv")
  write_mat(fit$stdev.unscaled, dir, "lmfit_stdev_unscaled.csv")
  write_mat(cbind(sigma = fit$sigma, df_residual = fit$df.residual, Amean = fit$Amean),
            dir, "lmfit_scalars.csv")
  write_mat(fit$cov.coefficients, dir, "lmfit_cov_coefficients.csv")
  cat(sprintf("%s: %d genes x %d arrays, %d blocks, consensus %.6f, %d NA rho\n", id,
              nrow(exprs), ncol(exprs), length(unique(block)), dc$consensus.correlation,
              sum(is.na(dc$atanh.correlations))))
}

# ---- techrep_synth / techrep_synth_na ------------------------------------------------------
set.seed(20261003)
ngenes <- 400
cond <- factor(rep(c("A", "B", "C", "D"), each = 6))
bio <- sprintf("S%02d", rep(1:12, each = 2))
mu <- rnorm(ngenes, 8, 1.5)
effect <- matrix(0, ngenes, 4)
de <- sample(ngenes, 80)
effect[de, 2:4] <- rnorm(80 * 3, 0, 1)
sd_bio <- sqrt(rgamma(ngenes, 4, 40))   # biological sample variance
sd_tech <- sqrt(rgamma(ngenes, 4, 80))  # technical replicate variance
bio_eff <- matrix(rnorm(ngenes * 12), ngenes, 12) * sd_bio
exprs <- mu + effect[, as.integer(cond)] + bio_eff[, rep(1:12, each = 2)] +
  matrix(rnorm(ngenes * 24), ngenes, 24) * sd_tech
colnames(exprs) <- sprintf("%s_%s_t%d", cond, bio, rep(1:2, 12))

design0 <- model.matrix(~0 + cond)
write_case("techrep_synth", exprs, design0, bio)

na_exprs <- exprs
na_exprs[matrix(runif(length(exprs)) < 0.12, ngenes)] <- NA
na_exprs[1, ] <- NA                         # all missing
na_exprs[2, -(1:4)] <- NA                   # two blocks, too few observations to fit
na_exprs[3, -c(1, 3, 5, 7, 9, 11)] <- NA    # one observation per block
na_exprs[4, -(1:2)] <- NA                   # one block only
na_exprs[5, c(TRUE, FALSE)] <- NA           # every block of size one
design1 <- model.matrix(~cond)
write_case("techrep_synth_na", na_exprs, design1, bio)

# ---- techrep_mixed -------------------------------------------------------------------------
set.seed(20261004)
ngenes <- 300
reps <- rep(c(3, 2, 3, 2), 4)
bio_id <- rep(seq_along(reps), reps)
cond <- factor(c("A", "B", "C", "D")[(bio_id - 1) %/% 4 + 1])
bio <- sprintf("S%02d", bio_id)
n <- length(bio)
mu <- rnorm(ngenes, 8, 1.5)
effect <- matrix(0, ngenes, 4)
de <- sample(ngenes, 60)
effect[de, 2:4] <- rnorm(60 * 3)
sd_bio <- sqrt(rgamma(ngenes, 4, 40))
sd_tech <- sqrt(rgamma(ngenes, 4, 80))
bio_eff <- matrix(rnorm(ngenes * 16), ngenes, 16) * sd_bio
exprs <- mu + effect[, as.integer(cond)] + bio_eff[, bio_id] +
  matrix(rnorm(ngenes * n), ngenes, n) * sd_tech
exprs[1, ] <- 8 + rnorm(16)[bio_id] + rnorm(n, 0, 1e-4)
third <- unlist(lapply(reps, function(k) seq_len(k) == 3))
dev <- unlist(lapply(reps, function(k) { e <- rnorm(1); c(e, -e, 0)[seq_len(k)] }))
exprs[2, ] <- 8 + as.integer(cond) + dev + rnorm(n, 0, 1e-3)
exprs[2, third] <- NA
colnames(exprs) <- sprintf("%s_%s_t%d", cond, bio, sequence(reps))
write_case("techrep_mixed", exprs, model.matrix(~0 + cond), bio)

# ---- bojkova_pairs -------------------------------------------------------------------------
mdir <- "tests/golden/corpus/matrix/conditions_only"
rd <- function(f) {
  d <- read.csv(file.path(mdir, f), check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(d[, -1]); rownames(m) <- d$row; m
}
bx <- rd("input_log2.csv")
bd <- rd("design.csv")
bcond <- colnames(bd)[max.col(bd)]
bblock <- character(ncol(bx))
for (cl in unique(bcond)) {
  idx <- which(bcond == cl)
  bblock[idx] <- sprintf("%s_p%d", cl, (seq_along(idx) + 1) %/% 2)
}
write_case("bojkova_pairs", bx, bd, bblock)
