# Times R limma on one dataset written by bench_data.R and prints one tab-separated row per
# stage: impl, dataset, pipeline, stage, seconds (the median of `reps` runs, each after a gc).
# Pipelines: `load` reads the inputs only and is the memory baseline; `fit` is lmFit,
# contrasts.fit, eBayes, topTable for every contrast and decideTests; `block` is
# duplicateCorrelation and then lmFit with the block and its consensus correlation.
# examples/bench.rs is the Rust twin and does the same calls in the same order.
# Usage: Rscript bench.R <dataset dir> <pipeline> <reps>
suppressPackageStartupMessages(library(limma))
args <- commandArgs(TRUE)
dir <- args[1]
pipeline <- args[2]
reps <- as.integer(args[3])
dataset <- basename(dir)

read <- function(file) as.matrix(read.csv(file.path(dir, file), check.names = FALSE))
exprs <- read("input_log2.csv")
design <- read("design.csv")
contrasts <- read("contrasts.csv")
rownames(contrasts) <- colnames(design)
block <- read.csv(file.path(dir, "block.csv"))$block

stage <- function(name, f) {
  secs <- numeric(reps)
  for (i in seq_len(reps)) {
    gc()
    start <- Sys.time()
    out <- f()
    secs[i] <- as.numeric(Sys.time() - start, units = "secs")
  }
  cat(sprintf("R\t%s\t%s\t%s\t%.6f\n", dataset, pipeline, name, median(secs)))
  invisible(out)
}

if (pipeline == "fit") {
  fit <- stage("lmFit", function() lmFit(exprs, design))
  cfit <- stage("contrasts.fit", function() contrasts.fit(fit, contrasts))
  eb <- stage("eBayes", function() eBayes(cfit))
  stage("topTable", function() {
    for (j in seq_len(ncol(contrasts))) {
      topTable(eb, coef = j, number = Inf, sort.by = "B", confint = TRUE)
    }
  })
  stage("decideTests", function() decideTests(eb, method = "separate", p.value = 0.05))
} else if (pipeline == "block") {
  dc <- stage("duplicateCorrelation", function() duplicateCorrelation(exprs, design, block = block))
  stage("lmFit(block)", function() {
    lmFit(exprs, design, block = block, correlation = dc$consensus.correlation)
  })
}
