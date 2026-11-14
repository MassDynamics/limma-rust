# Writes the synthetic datasets scripts/benchmark times R and Rust on, one folder per dataset
# under the directory given first, each with input_log2.csv, design.csv, contrasts.csv and
# block.csv. Every size comes complete and with 10% of values missing at random, because
# missing values send lmFit down a per-protein loop. Fixed seeds, so the data is the same on
# every run.
# Usage: Rscript bench_data.R <out dir> <protein counts...>
args <- commandArgs(TRUE)
out <- args[1]
sizes <- as.integer(args[-1])

# 24 samples: four conditions of six, from 12 donors with two technical replicates each.
group <- factor(rep(c("A", "B", "C", "D"), each = 6))
block <- rep(1:12, each = 2)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
contrasts <- limma::makeContrasts(B - A, C - A, D - A, C - B, D - B, D - C, levels = design)

for (n in sizes) {
  for (missing in c(0, 0.1)) {
    set.seed(n + 1000 * missing)
    effect <- matrix(0, n, 4)
    hit <- sample(n, n %/% 10)
    effect[hit, ] <- rnorm(length(hit) * 4)
    donor <- matrix(rnorm(n * 12, sd = 0.5), n)[, block]
    y <- 20 + rnorm(n, sd = 2) + effect[, as.integer(group)] + donor +
      matrix(rnorm(n * 24, sd = 0.3), n)
    if (missing > 0) y[sample(length(y), round(missing * length(y)))] <- NA
    colnames(y) <- paste0("s", 1:24)

    dir <- file.path(out, sprintf("%06d_%s", n, if (missing > 0) "na10" else "complete"))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(y, file.path(dir, "input_log2.csv"), row.names = FALSE)
    write.csv(design, file.path(dir, "design.csv"), row.names = FALSE)
    write.csv(contrasts, file.path(dir, "contrasts.csv"), row.names = FALSE)
    write.csv(data.frame(block = block), file.path(dir, "block.csv"), row.names = FALSE)
  }
}
