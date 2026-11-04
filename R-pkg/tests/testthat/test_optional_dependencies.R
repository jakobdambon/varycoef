test_that("optional dependency checks report packages needed by a feature", {
  unavailable <- function(package) FALSE

  expect_error(
    varycoef:::require_optional_dependencies(
      c("glmnet", "optimParallel"),
      "parallel selection",
      is_installed = unavailable
    ),
    "parallel selection requires these optional packages: glmnet, optimParallel"
  )
})

test_that("optional dependency checks pass when dependencies are available", {
  available <- function(package) TRUE

  expect_no_error(
    varycoef:::require_optional_dependencies(
      c("glmnet", "optimParallel"),
      "parallel selection",
      is_installed = available
    )
  )
})
