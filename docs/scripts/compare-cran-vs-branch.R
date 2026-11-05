#!/usr/bin/env Rscript

run_worker <- function(args) {
  if (length(args) != 8L) {
    stop("Internal worker received invalid arguments.", call. = FALSE)
  }

  lib <- args[[1L]]
  n <- as.integer(args[[2L]])
  repetitions <- as.integer(args[[3L]])
  cores <- as.integer(args[[4L]])
  dataset_file <- args[[5L]]
  result_file <- args[[6L]]
  tolerance <- as.numeric(args[[7L]])
  tapering <- if (identical(args[[8L]], "none")) {
    NULL
  } else {
    as.numeric(args[[8L]])
  }

  .libPaths(c(lib, .libPaths()))
  library(varycoef)
  datasets <- readRDS(dataset_file)
  dataset <- datasets[[as.character(n)]]
  train <- dataset$train

  cluster <- NULL
  if (cores > 1L) {
    if (!requireNamespace("optimParallel", quietly = TRUE)) {
      stop(
        "Parallel runs require the optional package 'optimParallel'. ",
        "Install it or rerun with --cores 1.",
        call. = FALSE
      )
    }

    Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep))
    cluster <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
  }

  control <- SVC_mle_control(
    profileLik = TRUE,
    # CRAN 0.3.6's tapered multi-SVC path passes an extra parameter to
    # covariance functions; use "exp", which accepts that legacy call.
    cov.name = "exp",
    tapering = tapering,
    # CRAN 0.3.6 reads this legacy field instead of `tapering`.
    taper = tapering,
    parallel = if (is.null(cluster)) NULL else list(cl = cluster)
  )

  fit_once <- function() {
    SVC_mle(
      y = train$y,
      X = train$X,
      locs = train$locs,
      W = train$W,
      control = control
    )
  }

  elapsed <- numeric(repetitions)
  fit <- NULL
  for (i in seq_len(repetitions)) {
    gc()
    elapsed[[i]] <- unname(system.time(fit <- fit_once())["elapsed"])
  }

  prediction <- predict(
    fit,
    newlocs = train$locs,
    newX = train$X,
    newW = train$W
  )
  fixed_effects <- coef(fit)
  estimated_svc <- as.matrix(prediction[, paste0("SVC_", seq_along(fixed_effects)),
                                         drop = FALSE])
  estimated_svc <- sweep(estimated_svc, 2L, fixed_effects, FUN = "+")

  result <- list(
    version = as.character(utils::packageVersion("varycoef")),
    n = n,
    tapering = if (is.null(tapering)) "none" else tapering,
    times = elapsed,
    outputs = list(
      coefficients = unname(coef(fit)),
      covariance_parameters = unname(cov_par(fit)),
      fitted_values = unname(fitted(fit)$y.pred),
      residuals = unname(residuals(fit)),
      predicted_values = unname(prediction$y.pred)
    ),
    plot_data = list(
      locs = drop(train$locs),
      true_svc = train$beta,
      estimated_svc = estimated_svc,
      response = train$y,
      fitted_response = prediction$y.pred
    ),
    tolerance = tolerance
  )
  saveRDS(result, result_file)
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) {
    stop("Could not determine the script path.", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

usage <- function() {
  cat(
    "Compare the CRAN varycoef release with the package in this checkout.\n\n",
    "Usage:\n",
    "  Rscript docs/scripts/compare-cran-vs-branch.R [options]\n\n",
    "Options:\n",
    "  --repetitions N  Timed fits per version and size (default: 10; minimum: 2)\n",
    "  --sizes N,N,N    Three increasing training sizes (default: 200,500,1000)\n",
    "  --tapering MODE  both, none, or a numeric taper range (default: both)\n",
    "                   The default compares untapered and tapered fits.\n",
    "  --tolerance X    Scaled numerical tolerance (default: 2e-3)\n",
    "  --cores N       Workers for optimParallel; 1 uses serial optim (default: 1)\n",
    "  --repo URL       CRAN mirror (default: https://cloud.r-project.org)\n",
    "  --output FILE   CSV report path (default: docs/benchmark-results.csv)\n",
    "  --plot-dir DIR  Directory for comparison PNGs (default: docs/images/benchmark)\n",
    "  --help           Show this help\n",
    sep = ""
  )
}

parse_options <- function(args) {
  options <- list(
    repetitions = 10L,
    sizes = c(200L, 500L, 1000L),
    taperings = c(NA_real_, 10),
    tolerance = 2e-3,
    cores = 1L,
    repo = "https://cloud.r-project.org",
    output = "docs/benchmark-results.csv",
    plot_dir = "docs/images/benchmark"
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (identical(key, "--help")) {
      usage()
      quit(save = "no", status = 0L)
    }
    if (!key %in% c(
      "--repetitions", "--sizes", "--tapering", "--tolerance", "--cores", "--repo",
      "--output", "--plot-dir"
    ) || i == length(args)) {
      usage()
      stop("Invalid or incomplete command-line option: ", key,
           call. = FALSE)
    }

    value <- args[[i + 1L]]
    if (identical(key, "--repetitions")) {
      options$repetitions <- suppressWarnings(as.integer(value))
    } else if (identical(key, "--sizes")) {
      options$sizes <- suppressWarnings(as.integer(strsplit(value, ",",
                                                            fixed = TRUE)[[1L]]))
    } else if (identical(key, "--tapering")) {
      options$taperings <- if (identical(value, "both")) {
        c(NA_real_, 10)
      } else if (identical(value, "none")) {
        NA_real_
      } else {
        numeric_value <- suppressWarnings(as.numeric(value))
        if (length(numeric_value) != 1L || is.na(numeric_value)) {
          stop("'--tapering' must be 'both', 'none', or a positive finite range.",
               call. = FALSE)
        }
        numeric_value
      }
    } else if (identical(key, "--tolerance")) {
      options$tolerance <- suppressWarnings(as.numeric(value))
    } else if (identical(key, "--cores")) {
      options$cores <- suppressWarnings(as.integer(value))
    } else if (identical(key, "--repo")) {
      options$repo <- value
    } else if (identical(key, "--output")) {
      options$output <- value
    } else {
      options$plot_dir <- value
    }
    i <- i + 2L
  }

  if (is.na(options$repetitions) || options$repetitions < 2L) {
    stop("'--repetitions' must be an integer of at least 2.", call. = FALSE)
  }
  if (length(options$sizes) != 3L || anyNA(options$sizes) ||
      any(options$sizes < 10L) || is.unsorted(options$sizes, strictly = TRUE)) {
    stop("'--sizes' must contain three increasing integers, each at least 10.",
         call. = FALSE)
  }
  if (!is.finite(options$tolerance) || options$tolerance <= 0) {
    stop("'--tolerance' must be a positive finite number.", call. = FALSE)
  }
  if (length(options$taperings) == 0L ||
      any(!is.na(options$taperings) &
          (!is.finite(options$taperings) | options$taperings <= 0))) {
    stop("'--tapering' must be 'both', 'none', or a positive finite range.",
         call. = FALSE)
  }
  if (is.na(options$cores) || options$cores < 1L) {
    stop("'--cores' must be a positive integer.", call. = FALSE)
  }
  options
}

run_subprocess <- function(rscript, args) {
  output <- system2(
    rscript,
    args = vapply(args, shQuote, character(1)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (status != 0L) {
    cat(paste(output, collapse = "\n"), "\n", file = stderr())
    stop("Subprocess failed with exit status ", status, ".", call. = FALSE)
  }
  invisible(output)
}

compare_outputs <- function(cran, branch, tolerance) {
  metrics <- names(cran$outputs)
  comparisons <- lapply(metrics, function(metric) {
    expected <- cran$outputs[[metric]]
    actual <- branch$outputs[[metric]]
    if (length(expected) != length(actual)) {
      return(list(matches = FALSE, max_scaled_difference = Inf))
    }
    difference <- if (length(expected)) {
      scaled_difference <- abs(expected - actual) / pmax(1, abs(expected))
      if (any(!is.finite(scaled_difference))) {
        Inf
      } else {
        max(scaled_difference)
      }
    } else {
      0
    }
    list(
      matches = is.finite(difference) && difference <= tolerance,
      max_scaled_difference = difference
    )
  })
  names(comparisons) <- metrics
  comparisons
}

write_comparison_plot <- function(n, runs, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(plot_dir, sprintf("svc-comparison-n-%d.png", n))
  scenarios <- unlist(lapply(names(runs), function(mode) {
    c(paste0("CRAN / ", mode), paste0("Branch / ", mode))
  }))
  grDevices::png(file, width = 1800, height = 450 * length(scenarios), res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(
    mfrow = c(length(scenarios), 3L),
    mar = c(3.5, 3.8, 2.5, 1),
    oma = c(1, 1, 3, 1)
  )
  on.exit(graphics::par(old_par), add = TRUE)
  row <- 0L
  for (mode in names(runs)) {
    for (version in c("cran", "branch")) {
      row <- row + 1L
      result <- runs[[mode]][[version]]
      plot_data <- result$plot_data
      location <- plot_data$locs
      row_title <- scenarios[[row]]

      for (j in seq_len(ncol(plot_data$true_svc))) {
        true <- plot_data$true_svc[, j]
        estimate <- plot_data$estimated_svc[, j]
        graphics::plot(
          location, true,
          type = "l",
          lwd = 2,
          xlab = "Location",
          ylab = sprintf("SVC %d", j),
          main = row_title,
          ylim = range(c(true, estimate))
        )
        graphics::lines(location, estimate, col = "firebrick", lwd = 1.5,
                        lty = 2)
        graphics::legend(
          "topright",
          legend = c("True", "Estimated"),
          col = c("black", "firebrick"),
          lty = c(1, 2),
          lwd = c(2, 1.5),
          bty = "n",
          cex = 0.8
        )
      }

      graphics::plot(
        location, plot_data$response,
        pch = 16,
        cex = 0.55,
        xlab = "Location",
        ylab = "Response",
        main = row_title
      )
      graphics::lines(
        location, plot_data$fitted_response,
        col = "firebrick",
        lwd = 1.5,
        lty = 2
      )
      graphics::legend(
        "topright",
        legend = c("True response", "Fitted response"),
        col = c("black", "firebrick"),
        pch = c(16, NA),
        lty = c(NA, 2),
        lwd = c(NA, 1.5),
        bty = "n",
        cex = 0.8
      )
    }
  }
  graphics::mtext(
    sprintf(
      "True and estimated SVCs and response (training observations, n = %d)",
      n
    ),
    outer = TRUE,
    side = 3,
    line = 1
  )
  invisible(file)
}

main <- function(args) {
  options <- parse_options(args)
  script <- script_path()
  package_dir <- normalizePath(file.path(dirname(script), "..", "..", "R-pkg"),
                               mustWork = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  r_command <- file.path(R.home("bin"), "R")

  if (options$cores > 1L && !requireNamespace("optimParallel", quietly = TRUE)) {
    stop(
      "Parallel benchmarking requires 'optimParallel' in the current R library. ",
      "Install it or use --cores 1.",
      call. = FALSE
    )
  }

  temp_root <- tempfile("varycoef-comparison-")
  dir.create(temp_root)
  on.exit(unlink(temp_root, recursive = TRUE), add = TRUE)
  cran_lib <- file.path(temp_root, "cran-library")
  branch_lib <- file.path(temp_root, "branch-library")
  dir.create(cran_lib)
  dir.create(branch_lib)

  .libPaths(c(cran_lib, branch_lib, .libPaths()))
  cran_version <- "0.3.6"
  if (!requireNamespace("remotes", quietly = TRUE)) {
    utils::install.packages(
      "remotes",
      lib = cran_lib,
      repos = options$repo
    )
  }
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("Could not install the 'remotes' package needed for the CRAN archive.",
         call. = FALSE)
  }

  message(sprintf(
    "Installing the pinned CRAN varycoef version %s...",
    cran_version
  ))
  remotes::install_version(
    package = "varycoef",
    version = cran_version,
    lib = cran_lib,
    repos = options$repo,
    type = "source",
    dependencies = NA,
    upgrade = "never"
  )
  if (!dir.exists(file.path(cran_lib, "varycoef"))) {
    stop("The archived CRAN varycoef package was not installed.",
         call. = FALSE)
  }
  installed_cran_version <- as.character(
    utils::packageVersion("varycoef", lib.loc = cran_lib)
  )
  if (!identical(installed_cran_version, cran_version)) {
    stop(
      sprintf(
        "Expected CRAN varycoef %s, but installed %s.",
        cran_version,
        installed_cran_version
      ),
      call. = FALSE
    )
  }

  message("Installing the package from this checkout...")
  install_output <- system2(
    r_command,
    args = c(
      "CMD", "INSTALL", "-l", shQuote(branch_lib),
      shQuote(package_dir)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  install_status <- attr(install_output, "status")
  if (!is.null(install_status) && install_status != 0L) {
    cat(paste(install_output, collapse = "\n"), "\n", file = stderr())
    stop("Could not install this checkout into the temporary library.",
         call. = FALSE)
  }
  if (!dir.exists(file.path(branch_lib, "varycoef"))) {
    stop("The checkout's varycoef package was not installed.", call. = FALSE)
  }

  branch_namespace <- loadNamespace("varycoef", lib.loc = branch_lib)
  sample_svc_data <- get("sample_SVCdata", envir = branch_namespace)
  maximum_n <- max(options$sizes)
  set.seed(20261105L)
  synthetic_data <- sample_svc_data(
    df.pars = data.frame(
      var = c(2, 1),
      scale = c(3, 1),
      mean = c(1, 2)
    ),
    nugget.sd = 0.5,
    locs = seq_len(maximum_n + 5L),
    cov.name = "exp"
  )
  datasets <- lapply(seq_along(options$sizes), function(i) {
    n <- options$sizes[[i]]
    train_id <- seq_len(n)
    list(
      train = list(
        y = synthetic_data$y[train_id],
        X = synthetic_data$X[train_id, , drop = FALSE],
        W = synthetic_data$X[train_id, , drop = FALSE],
        locs = synthetic_data$locs[train_id, , drop = FALSE],
        beta = synthetic_data$beta[train_id, , drop = FALSE]
      )
    )
  })
  names(datasets) <- as.character(options$sizes)
  dataset_file <- file.path(temp_root, "datasets.rds")
  saveRDS(datasets, dataset_file)
  rm(branch_namespace, sample_svc_data, synthetic_data, datasets)
  gc()

  results <- list()
  for (i in seq_along(options$sizes)) {
    n <- options$sizes[[i]]
    runs_by_tapering <- list()
    for (tapering in options$taperings) {
      taper_label <- if (is.na(tapering)) "none" else as.character(tapering)
      mode <- if (is.na(tapering)) "none" else paste0("range-", taper_label)
      message(sprintf(
        "Benchmarking n=%d, tapering=%s (%d timed fits per version)...",
        n, taper_label, options$repetitions
      ))
      order <- if (i %% 2L) c("cran", "branch") else c("branch", "cran")
      runs <- list()
      for (version in order) {
        lib <- if (identical(version, "cran")) cran_lib else branch_lib
        result_file <- file.path(
          temp_root, sprintf("%s-%d-%s.rds", version, n, mode)
        )
        taper_arg <- if (is.na(tapering)) "none" else as.character(tapering)
        run_subprocess(
          rscript,
          c(
            script, "--worker", lib, n, options$repetitions, options$cores,
            dataset_file, result_file, options$tolerance, taper_arg
          )
        )
        runs[[version]] <- readRDS(result_file)
      }

      comparisons <- compare_outputs(
        runs$cran, runs$branch, options$tolerance
      )
      cran_median <- median(runs$cran$times)
      branch_median <- median(runs$branch$times)
      row <- data.frame(
        n = n,
        cran_version = runs$cran$version,
        branch_version = runs$branch$version,
        repetitions = options$repetitions,
        optimizer_workers = options$cores,
        tapering_range = if (is.na(tapering)) "none" else taper_label,
        tolerance = options$tolerance,
        cran_median_seconds = cran_median,
        branch_median_seconds = branch_median,
        speedup_percent = 100 * (cran_median - branch_median) / cran_median,
        results_match = all(vapply(comparisons, `[[`, logical(1), "matches")),
        max_scaled_difference = max(vapply(
          comparisons, `[[`, numeric(1), "max_scaled_difference"
        ))
      )
      results[[length(results) + 1L]] <- row
      runs_by_tapering[[taper_label]] <- runs
      message(sprintf(
        paste0(
          "  numerical match: %s; CRAN median %.3f s; branch median %.3f s; ",
          "speedup %.1f%%"
        ),
        row$results_match, cran_median, branch_median, row$speedup_percent
      ))
    }

    plot_file <- write_comparison_plot(n, runs_by_tapering, options$plot_dir)
    message("  plot: ", normalizePath(plot_file, mustWork = FALSE))
  }

  report <- do.call(rbind, results)
  utils::write.csv(report, options$output, row.names = FALSE)
  cat("\nComparison summary:\n")
  print(report, row.names = FALSE)
  cat("\nCSV report: ", normalizePath(options$output, mustWork = FALSE), "\n",
      sep = "")
  if (!all(report$results_match)) {
    stop("At least one result comparison exceeded the requested tolerance.",
         call. = FALSE)
  }
  if (all(report$speedup_percent > 0)) {
    cat("The branch was faster at all tested dataset sizes.\n")
  } else {
    cat(
      "The branch was not faster at every size; consider more repetitions ",
      "or larger sizes to reduce timing noise.\n",
      sep = ""
    )
  }
}

trailing <- commandArgs(trailingOnly = TRUE)
if (length(trailing) && identical(trailing[[1L]], "--worker")) {
  run_worker(trailing[-1L])
} else {
  main(trailing)
}
