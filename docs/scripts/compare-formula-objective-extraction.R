#!/usr/bin/env Rscript

# Run from the repository root:
# Rscript docs/scripts/compare-formula-objective-extraction.R
# The baseline is the master commit before this fix. Git, R, and the package's
# required dependencies must be available; no packages are downloaded.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) && identical(args[[1L]], "--worker")) {
  if (length(args) != 3L) stop("Invalid worker arguments.", call. = FALSE)
  .libPaths(c(args[[2L]], .libPaths()))
  library(varycoef)
  data("SVCdata", package = "varycoef")
  id <- seq_len(12L)
  df <- data.frame(y = SVCdata$y[id], x = SVCdata$X[id, 2L])
  X <- model.matrix(y ~ x, data = df)
  W <- model.matrix(~ 1, data = df)
  locs <- SVCdata$locs[id]
  control <- SVC_mle_control(extract_fun = TRUE, profileLik = TRUE)
  pars <- c(0.5, 0.5, 0.5)

  matrix_obj <- SVC_mle(df$y, X, locs, W = W, control = control)
  matrix_value <- do.call(matrix_obj$obj_fun, c(list(x = pars), matrix_obj$args))
  formula_result <- tryCatch(
    {
      obj <- SVC_mle(
        y ~ x, data = df, RE_formula = ~ 1, locs = locs, control = control
      )
      if (!inherits(obj, "SVC_obj_fun")) {
        stop("Formula call did not return an SVC_obj_fun object.")
      }
      list(
        status = "success",
        value = do.call(obj$obj_fun, c(list(x = pars), obj$args))
      )
    },
    error = function(e) list(status = "error", message = conditionMessage(e))
  )
  saveRDS(
    list(
      version = as.character(utils::packageVersion("varycoef")),
      matrix_value = matrix_value,
      formula = formula_result
    ),
    args[[3L]]
  )
} else {
  if (length(args)) stop("This script does not accept arguments.", call. = FALSE)

  main <- function() {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (!length(file_arg)) stop("Could not determine the script path.", call. = FALSE)
    script <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
    root <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
    baseline <- "cd3930bce5813a7a54bc2fe8eeab4f167cb3026f"
    work <- tempfile("varycoef-formula-comparison-")
    dir.create(work)
    on.exit(unlink(work, recursive = TRUE), add = TRUE)
    dir.create(file.path(work, "baseline"))
    dir.create(file.path(work, "baseline-lib"))
    dir.create(file.path(work, "current-lib"))

    run <- function(command, arguments) {
      output <- system2(
        command, args = vapply(arguments, shQuote, character(1)),
        stdout = TRUE, stderr = TRUE
      )
      status <- attr(output, "status")
      if (!is.null(status) && status != 0L) {
        stop(paste(c(paste(command, "failed:"), output), collapse = "\n"),
             call. = FALSE)
      }
      invisible(output)
    }

    archive <- file.path(work, "baseline.tar")
    run("git", c("-C", root, "archive", "--format=tar",
                 paste0("--output=", archive), baseline, "R-pkg"))
    utils::untar(archive, exdir = file.path(work, "baseline"))

    r <- file.path(R.home("bin"), "R")
    rscript <- file.path(R.home("bin"), "Rscript")
    sources <- c(file.path(work, "baseline", "R-pkg"), file.path(root, "R-pkg"))
    libs <- file.path(work, c("baseline-lib", "current-lib"))
    outputs <- file.path(work, c("baseline.rds", "current.rds"))
    for (i in seq_along(sources)) {
      run(r, c("CMD", "INSTALL", "--library", libs[[i]], sources[[i]]))
      run(rscript, c(script, "--worker", libs[[i]], outputs[[i]]))
    }
    before <- readRDS(outputs[[1L]])
    after <- readRDS(outputs[[2L]])
    if (before$formula$status != "error") {
      stop("The recorded baseline unexpectedly accepted the formula call.",
           call. = FALSE)
    }
    if (after$formula$status != "success") {
      stop("The current formula call failed: ", after$formula$message,
           call. = FALSE)
    }
    if (!isTRUE(all.equal(after$formula$value, after$matrix_value,
                          tolerance = 1e-10))) {
      stop("The current formula and matrix objectives differ.", call. = FALSE)
    }
    if (!isTRUE(all.equal(before$matrix_value, after$matrix_value,
                          tolerance = 1e-10))) {
      stop("The matrix objective differs between versions.", call. = FALSE)
    }

    cat("Baseline master commit: ", substr(baseline, 1L, 7L),
        " (version ", before$version, ")\n",
        "Current checkout version: ", after$version, "\n",
        "Baseline matrix objective: ", before$matrix_value, "\n",
        "Baseline formula call: error: ", before$formula$message, "\n",
        "Current matrix objective: ", after$matrix_value, "\n",
        "Current formula objective: ", after$formula$value, "\n",
        "Current formula and matrix objectives agree.\n", sep = "")
  }
  main()
}
