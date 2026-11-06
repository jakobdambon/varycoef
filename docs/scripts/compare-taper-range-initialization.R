#!/usr/bin/env Rscript

# Compare the tapered range initializer on master before the fix with the
# current implementation. Run from any directory with:
# Rscript docs/scripts/compare-taper-range-initialization.R

script_arg <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (!length(script_arg)) {
  stop("Could not determine the script path.", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1L]]))
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."))

master_commit <- "c8ff4196768ffcfa787e282a1f3e9e7b44253829"
master_source <- system2(
  "git",
  args = c(
    "-C", shQuote(repo_root),
    "show", shQuote(paste0(master_commit, ":R-pkg/R/SVC_mle.R"))
  ),
  stdout = TRUE,
  stderr = TRUE
)
if (!is.null(attr(master_source, "status")) &&
    attr(master_source, "status") != 0L) {
  stop(
    "Could not read the recorded master version. Fetch the repository history ",
    "and try again.",
    call. = FALSE
  )
}
if (!any(grepl("median(lower.tri(d@entries))", master_source, fixed = TRUE))) {
  stop("The recorded master source does not contain the expected old initializer.",
       call. = FALSE)
}

current_source <- readLines(
  file.path(repo_root, "R-pkg/R/SVC_mle.R"),
  warn = FALSE
)
if (!any(grepl("median_sparse_distance(d)", current_source, fixed = TRUE))) {
  stop("The current SVC_mle source does not use the fixed initializer.",
       call. = FALSE)
}

helper_env <- new.env(parent = baseenv())
source(file.path(repo_root, "R-pkg/R/utils.R"), local = helper_env)
if (!exists("median_sparse_distance", envir = helper_env, inherits = FALSE)) {
  stop("The current sparse-distance helper could not be loaded.",
       call. = FALSE)
}

current_commit <- system2(
  "git",
  args = c("-C", shQuote(repo_root), "rev-parse", "--short", "HEAD"),
  stdout = TRUE,
  stderr = TRUE
)
cat(
  "Master reference: ", substr(master_commit, 1L, 7L),
  "\nCurrent checkout: ", current_commit[[1L]],
  "\n\n",
  sep = ""
)
cat(
  "Example: six equally spaced locations; taper range scales with the ",
  "coordinate units.\n",
  "The default range start/lower/upper are median distance / 4, ",
  "median distance / 1000, and 10 * median distance.\n\n",
  sep = ""
)

scales <- c(1, 10)
results <- lapply(scales, function(scale) {
  locs <- matrix(scale * seq_len(6L), ncol = 1L)
  distances <- spam::nearest.dist(locs, delta = 3 * scale)

  master_median <- stats::median(base::lower.tri(distances@entries))
  current_median <- helper_env$median_sparse_distance(distances)

  data.frame(
    coordinate_scale = scale,
    version = c("master (before fix)", "current (fixed)"),
    median_distance = c(master_median, current_median),
    range_start = c(master_median, current_median) / 4,
    lower_bound = c(master_median, current_median) / 1000,
    upper_bound = 10 * c(master_median, current_median)
  )
})
comparison <- do.call(rbind, results)
rownames(comparison) <- NULL

expected_medians <- c(1, 2, 1, 20)
if (!isTRUE(all.equal(comparison$median_distance, expected_medians))) {
  stop("The comparison did not produce the expected scale-dependent medians.",
       call. = FALSE)
}
print(comparison, digits = 5L, row.names = FALSE)

cat(
  "\nOn master, lower.tri() was applied to the sparse distance entries ",
  "vector, producing a logical mask whose median is always 1 here. ",
  "Consequently, the default range start and bounds do not change when ",
  "coordinates are rescaled. The fixed implementation takes the median ",
  "of retained positive pairwise distances, so its defaults scale with ",
  "the coordinates.\n",
  sep = ""
)
