test_that("prediction functions works with formula and matrix calls", {
  ## ---- toy example (from SVC_mle help file) ----
  # We use the sampled, i.e., onde dimensional SVCs
  data(SVCdata)
  # sub-sample data to have feasible run time for example
  set.seed(123)
  id <- sample(length(SVCdata$locs), 50)
  
  ## SVC_mle call with matrix arguments
  fit_mat <- with(SVCdata, SVC_mle(
    y[id], X[id, ], locs[id], 
    control = SVC_mle_control(profileLik = TRUE, cov.name = "mat32")))
  
  ## SVC_mle call with formula
  df <- with(SVCdata, data.frame(y = y[id], X = X[id, -1]))
  fit_form <- SVC_mle(
    y ~ X, data = df, locs = SVCdata$locs[id], 
    control = SVC_mle_control(profileLik = TRUE, cov.name = "mat32")
  )
  
  ## ---- predictions ----
  newdata <- data.frame(X = 3:4)
  newlocs <- 1:2
  newX <- matrix(c(1, 1, 3:4), ncol = 2)
  # only predicting SVC and response 
  pred_mat <- predict(fit_mat, newX = newX, newW = newX, newlocs = newlocs)
  pred_form <- predict(fit_form, newdata = newdata, newlocs = newlocs)
  
  expect_equal(pred_mat, pred_form, tolerance = 1e-10)
  
  # only predicting SVCs 
  pred_mat <- predict(fit_mat, newlocs = newlocs)
  pred_form <- predict(fit_form, newlocs = newlocs)
  
  expect_equal(pred_mat, pred_form, tolerance = 1e-10)
  # check warning for overwriting arguments
  expect_warning(
    predict(fit_form, newdata = newdata, newX = newX, newlocs = newlocs))
  expect_warning(
    predict(fit_form, newdata = newdata, newW = newX, newlocs = newlocs))
})