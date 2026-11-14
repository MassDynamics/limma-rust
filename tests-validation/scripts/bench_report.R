# Turns scripts/benchmark's benchmark.tsv into benchmark.md (summary and tables) and
# benchmark.svg (heatmaps of how many times faster and leaner Rust is than R).
# Usage: Rscript bench_report.R <out dir> <machine> <rust version> <reps>
args <- commandArgs(TRUE)
out <- args[1]
machine <- args[2]
rust <- args[3]
reps <- args[4]

x <- read.delim(file.path(out, "benchmark.tsv"), stringsAsFactors = FALSE)
datasets <- sort(unique(x$dataset))
proteins <- function(d) format(as.integer(sub("_.*", "", d)), big.mark = ",")
missing <- function(d) if (grepl("na10", d)) "10% missing" else "complete"
label <- function(d) paste(proteins(d), missing(d))
stages <- c("lmFit", "contrasts.fit", "eBayes", "topTable", "decideTests",
            "duplicateCorrelation", "lmFit(block)")
wide <- function(rows) unclass(xtabs(value ~ stage + dataset, rows))

times <- x[x$stage != "peak_rss_mb", ]
r_time <- wide(times[times$impl == "R", ])[stages, datasets]
rust_time <- wide(times[times$impl == "Rust", ])[stages, datasets]
speedup <- r_time / rust_time

mem <- x[x$stage == "peak_rss_mb", ]
peak <- function(impl, pipeline) {
  wide(transform(mem[mem$impl == impl & mem$pipeline == pipeline, ], stage = pipeline))[1, datasets]
}
mem_ratio <- rbind(
  "peak memory, fit" = peak("R", "fit") / peak("Rust", "fit"),
  "peak memory, block" = peak("R", "block") / peak("Rust", "block")
)

# The chart: one cell per stage and dataset, green where Rust wins and red where R does, on a
# log scale so 10x and 0.1x sit equally far from white.
times_x <- function(v) {
  ifelse(v >= 10, sprintf("%.0fx", v), ifelse(v >= 1, sprintf("%.1fx", v), sprintf("%.2fx", v)))
}
pal <- colorRampPalette(c("#d73027", "#f7f7f7", "#1a9850"))(101)
heat <- function(m, title, top) {
  shade <- pal[pmin(pmax(round((log10(m) + 2) / 4 * 100) + 1, 1), 101)]
  cell_w <- 110; cell_h <- 26; left <- 170
  svg <- c(sprintf('<text x="10" y="%d" font-weight="bold">%s</text>', top + 16, title))
  for (j in seq_along(datasets)) {
    svg <- c(svg, sprintf('<text x="%d" y="%d" text-anchor="middle">%s</text><text x="%d" y="%d" text-anchor="middle" fill="#555">%s</text>',
                          left + (j - 0.5) * cell_w, top + 40, proteins(datasets[j]),
                          left + (j - 0.5) * cell_w, top + 54, missing(datasets[j])))
  }
  for (i in seq_len(nrow(m))) {
    y <- top + 62 + (i - 1) * cell_h
    svg <- c(svg, sprintf('<text x="%d" y="%d" text-anchor="end">%s</text>', left - 8, y + 17, rownames(m)[i]))
    for (j in seq_len(ncol(m))) {
      svg <- c(svg, sprintf('<rect x="%d" y="%d" width="%d" height="%d" fill="%s" stroke="white" stroke-width="2"/><text x="%d" y="%d" text-anchor="middle">%s</text>',
                            left + (j - 1) * cell_w, y, cell_w, cell_h, shade[(j - 1) * nrow(m) + i],
                            left + (j - 0.5) * cell_w, y + 17, times_x(m[i, j])))
    }
  }
  list(svg = svg, bottom = top + 62 + nrow(m) * cell_h)
}
time_panel <- heat(speedup, "Time: R time / Rust time (green: Rust faster)", 0)
mem_panel <- heat(mem_ratio, "Memory: R peak / Rust peak", time_panel$bottom + 20)
writeLines(c(
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" font-family="sans-serif" font-size="12">',
          170 + 110 * length(datasets) + 20, mem_panel$bottom + 10),
  '<rect width="100%" height="100%" fill="white"/>', time_panel$svg, mem_panel$svg, "</svg>"
), file.path(out, "benchmark.svg"))

# The markdown.
fmt_s <- function(v) ifelse(v >= 1, sprintf("%.2f", v), ifelse(v >= 0.01, sprintf("%.3f", v), sprintf("%.1e", v)))
best <- which(speedup == max(speedup), arr.ind = TRUE)[1, ]
worst <- which(speedup == min(speedup), arr.ind = TRUE)[1, ]
lines <- c(
  "# limma-rust against R limma: speed and memory", "",
  sprintf("- Machine: %s", machine),
  sprintf("- R: %s, limma %s, %s", R.version.string, packageVersion("limma"),
          basename(sessionInfo()$BLAS)),
  sprintf("- Rust: %s, release build", rust),
  "- Data: synthetic, 24 samples (four conditions of six, 12 donors with two technical replicates each), complete and with 10% of values missing at random",
  sprintf("- Rust was faster in %d of %d stage and dataset pairs; the median speed-up is %s, from %s (%s, %s) to %s (%s, %s).",
          sum(speedup > 1), length(speedup), times_x(median(speedup)),
          times_x(speedup[worst[1], worst[2]]), stages[worst[1]], label(datasets[worst[2]]),
          times_x(speedup[best[1], best[2]]), stages[best[1]], label(datasets[best[2]])),
  sprintf("- Peak memory: R %.0f to %.0f MB, Rust %.0f to %.0f MB.",
          min(mem$value[mem$impl == "R" & mem$pipeline != "load"]),
          max(mem$value[mem$impl == "R" & mem$pipeline != "load"]),
          min(mem$value[mem$impl == "Rust" & mem$pipeline != "load"]),
          max(mem$value[mem$impl == "Rust" & mem$pipeline != "load"])),
  "", "![Rust against R, time and memory](benchmark.svg)", "",
  "## How to read this", "",
  "Both sides run the same limma calls on the same files, natively on the machine above and single-threaded. `fit` is lmFit, contrasts.fit, eBayes, topTable for each of the six contrasts and decideTests; `block` is duplicateCorrelation and then lmFit with the block and its consensus correlation.",
  "",
  sprintf("- **Time** is the median of %s runs of each stage inside one process, with the inputs already in memory, so start-up and reading files are left out. The speed-up is R's time divided by Rust's: 10x means Rust took a tenth of the time, under 1x means R was faster.", reps),
  "- **Peak memory** is the most resident memory the whole process held, from `/usr/bin/time`, one process per pipeline. It includes the runtime itself, which for R is the interpreter and limma, so it is what a container running the step would need. **Added** is that peak minus the peak of a process that only reads the inputs, which is closer to what the computation itself allocated.",
  "- Missing values matter because they send R's lmFit from one matrix solve for all proteins to a loop that fits each protein on its own.",
  "- These are one machine's numbers on one day and move by tens of percent between runs, so read the size of a ratio, not its last digit. Nothing fails on them.",
  "", "## Time (seconds)", "",
  "| dataset | stage | R | Rust | R / Rust |", "|---|---|---:|---:|---:|"
)
for (d in datasets) for (s in stages) {
  lines <- c(lines, sprintf("| %s | %s | %s | %s | %s |", label(d), s, fmt_s(r_time[s, d]),
                            fmt_s(rust_time[s, d]), times_x(speedup[s, d])))
}
lines <- c(lines, "", "## Memory (MB)", "",
           "| dataset | pipeline | R peak | Rust peak | R / Rust | R added | Rust added |",
           "|---|---|---:|---:|---:|---:|---:|")
for (d in datasets) for (p in c("fit", "block")) {
  rp <- peak("R", p)[d]; up <- peak("Rust", p)[d]
  lines <- c(lines, sprintf("| %s | %s | %.0f | %.1f | %s | %.0f | %.1f |", label(d), p, rp, up,
                            times_x(rp / up), rp - peak("R", "load")[d], up - peak("Rust", "load")[d]))
}
writeLines(lines, file.path(out, "benchmark.md"))
