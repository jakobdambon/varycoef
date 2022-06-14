test_that("SVC_mle call creates correct objects", {
  ## ---- toy example ----
  ## We use the sampled, i.e., onde dimensional SVCs
  data(SVCdata)
  # sub-sample data to have feasible run time for example
  set.seed(123)
  id <- sample(length(SVCdata$locs), 50)
  
  ## SVC_mle call with matrix arguments
  fit1 <- with(SVCdata, SVC_mle(
    y[id], X[id, ], locs[id], 
    control = SVC_mle_control(profileLik = TRUE, cov.name = "mat32")))
  
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