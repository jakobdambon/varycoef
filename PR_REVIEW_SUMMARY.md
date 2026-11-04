# PR review summary

## Review scope

This branch consolidates the package-improvement work on `review/package-quality`,
based on `master` at `ba8179a`. The implementation is a linear stack of seven
focused commits; the branch tip is `5679685` before this summary is added.

## Changes

- Added validation for model, selection, and prediction inputs, with clearer
  errors for invalid values and dimensions.
- Made selection-related dependencies optional, retaining selection for users
  who install those optional packages.
- Added an R CMD check workflow for R release and devel, including Pandoc setup
  for vignette checks.
- Replaced general covariance solves with triangular solves in model and
  selection paths.
- Added regression tests for validation, matrix/formula fit agreement,
  prediction consistency, predictive variance, and optional-dependency errors.
- Fixed predictive response covariance construction so it passes the covariance
  function the expected range and variance parameters.

## Dependency impact

The mandatory `Imports` list changes from 10 packages on `master` to 3 on this
branch (`methods`, `pbapply`, and `spam`). Seven packages used by selection or
parallel selection (`glmnet`, `lhs`, `mlr`, `mlrMBO`, `optimParallel`,
`ParamHelpers`, and `smoof`) move to `Suggests`. The selection feature remains
available when its relevant optional packages are installed; core MLE fitting
no longer requires them.

## Runtime measurement

Compared `master` and this branch in the same R 4.5.1/macOS ARM environment,
using all 500 observations in `SVCdata`, `SVC_mle()` with
`profileLik = TRUE` and `cov.name = "mat32"`, and three consecutive timed fits
per version:

| Version | Fit times (seconds) | Median |
|---|---:|---:|
| `master` | 132.104, 130.283, 130.611 | 130.611 s |
| `review/package-quality` | 100.432, 99.863, 100.625 | 100.432 s |

The branch median is 23.1% lower for this run. This is a local benchmark, not a
performance guarantee; runtime depends on hardware, R, and the input data.

## Evidence that core capabilities are retained

- A comparison of the loaded namespaces found the same 10 exported names and
  the same formal arguments for all exported functions on `master` and this
  branch.
- The complete `master` test suite passed: 48 assertions, 0 failures, 0
  warnings. The branch suite passed: 98 assertions, 0 failures, 0 warnings.
- Both versions fitted the same 500-row `SVCdata` input and returned 500
  predictions. Comparing the fits and predictions, the maximum scaled
  differences were approximately `0.00018` for coefficients, `0.00121` for
  covariance parameters, and `0.00017` for fitted and predicted values. All
  comparisons passed at tolerance `1e-3`; small numerical differences are
  expected from the changed linear-solve calculations.
- Package checking passed tests, examples, code, and documentation checks. The
  local check used `--no-build-vignettes` because Pandoc was unavailable, so it
  reported two vignette warnings; CI installs Pandoc and checks R release and
  devel.

The regression tests and namespace comparison show that the established public
API and core fit/predict workflows remain available. Selection dependencies
are optional rather than removed, and the new validation rejects invalid
inputs earlier.
