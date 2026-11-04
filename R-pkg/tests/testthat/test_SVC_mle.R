test_that("SVC_mle call creates correct objects", {
  ## ---- toy example ----
  ## We use the sampled, i.e., onde dimensional SVCs
  data(SVCdata)
  # sub-sample data to have feasible run time for example
  set.seed(123)
  id <- sample(length(SVCdata$locs), 50)
  
  ## SVC_mle call with matrix arguments
  
  # Test on the side: Internally, the spam.trivalues option is set to TRUE. 
  # We check that it is set back to the original value after the call.
  options(spam.trivalues = FALSE)
  expect_identical(FALSE, getOption("spam.trivalues"))
  
  fit1 <- with(SVCdata, SVC_mle(
    y[id], X[id, ], locs[id], 
    control = SVC_mle_control(profileLik = TRUE, cov.name = "mat32")))
  
  expect_identical(FALSE, getOption("spam.trivalues"))
  
  ## SVC_mle call with formula
  df <- with(SVCdata, data.frame(y = y[id], X = X[id, -1]))
  fit2 <- SVC_mle(
    y ~ X, data = df, locs = SVCdata$locs[id], 
    control = SVC_mle_control(profileLik = TRUE, cov.name = "mat32")
  )
  
  ## SVC_mle call with formula
  df <- with(SVCdata, data.frame(y = y[id], X = X[id, -1]))
  fit2_2 <- SVC_mle(
    y ~ X, data = df, RE_formula = ~ 1, locs = SVCdata$locs[id], 
    control = SVC_mle_control(profileLik = TRUE, cov.name = "mat32")
  )
  
  expect_null(fit1$formula, fit1$RE_formula)
  
  expect_identical(fit2$formula, y ~ X)
  expect_identical(fit2$RE_formula, y ~ X)
  
  expect_identical(fit2_2$formula, y ~ X)
  expect_identical(fit2_2$RE_formula, ~ 1)
})

test_that("SVC_mle validates matrix inputs before fitting", {
  y <- c(1, 2, 3)
  X <- cbind(1, c(0, 1, 2))
  locs <- c(0, 1, 2)

  expect_error(SVC_mle(y[-1], X, locs), "'X' must have 2 rows")
  expect_error(SVC_mle(y, X[, 1, drop = FALSE], locs,
                       W = matrix(1:2, ncol = 1)), "'W' must have 3 rows")
  expect_error(SVC_mle(y, X, locs[-1]), "'locs' must have 3 rows")
  expect_error(SVC_mle(c(1, NA, 3), X, locs), "'y' must not contain")
  expect_error(SVC_mle(y, cbind(X, Inf), locs), "'X' must not contain")
  expect_error(SVC_mle(y, X, c(0, NaN, 2)), "'locs' must not contain")
})

test_that("formula SVC_mle validates the response and row alignment", {
  data <- data.frame(y = c(1, 2, 3), x = 1:3)
  locs <- 1:3

  expect_error(
    SVC_mle(y ~ x, data = transform(data, y = c(1, NA, 3)), locs = locs),
    "'y' must not contain"
  )
  expect_error(
    SVC_mle(y ~ x, data = data, locs = locs[-1]),
    "'locs' must have 3 rows"
  )
  expect_error(
    SVC_mle(y ~ x, data = list(y = 1:3, x = 1:3), locs = locs),
    "'data' must be a data frame"
  )
})

test_that("matrix and formula fits agree on estimates and fitted values", {
  data(SVCdata)
  set.seed(123)
  id <- sample(length(SVCdata$locs), 50)
  y <- SVCdata$y[id]
  X <- SVCdata$X[id, , drop = FALSE]
  locs <- SVCdata$locs[id]
  df <- data.frame(y = y, X = X[, 2])
  control <- SVC_mle_control(
    profileLik = TRUE,
    cov.name = "mat32"
  )

  fit_matrix <- SVC_mle(y, X, locs, control = control)
  fit_formula <- SVC_mle(y ~ X, data = df, locs = locs, control = control)

  expect_equal(unname(coef(fit_matrix)), unname(coef(fit_formula)),
               tolerance = 1e-8)
  expect_equal(unname(cov_par(fit_matrix)), unname(cov_par(fit_formula)),
               tolerance = 1e-8)
  expect_equal(fitted(fit_matrix)$y.pred, fitted(fit_formula)$y.pred,
               tolerance = 1e-8)
  expect_equal(residuals(fit_matrix), residuals(fit_formula),
               tolerance = 1e-8)
})
