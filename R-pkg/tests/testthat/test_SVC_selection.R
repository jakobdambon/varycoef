selection_objective <- function() {
  data(SVCdata)
  id <- seq_len(12)
  SVC_mle(
    y = SVCdata$y[id],
    X = SVCdata$X[id, , drop = FALSE],
    locs = SVCdata$locs[id],
    W = SVCdata$X[id, 1, drop = FALSE],
    control = SVC_mle_control(
      extract_fun = TRUE,
      profileLik = TRUE,
      cov.name = "mat32"
    )
  )
}

test_that("SVC_selection_control accepts valid settings", {
  expect_type(SVC_selection_control(), "list")
  expect_identical(
    SVC_selection_control(method = "MBO", n.init = 3)$method,
    "MBO"
  )
  expect_identical(
    SVC_selection_control(
      r.lambda = c(1e-8, 1),
      CD.conv = list(N = 4L, delta = 1e-5, logLik = FALSE),
      adaptive = TRUE,
      hessian = TRUE
    )$CD.conv$logLik,
    FALSE
  )
})

test_that("SVC_selection_control rejects invalid settings", {
  expect_error(SVC_selection_control(r.lambda = c(1, 1)), "r.lambda")
  expect_error(SVC_selection_control(r.lambda = c(0, 1)), "r.lambda")
  expect_error(SVC_selection_control(n.lambda = 0), "n.lambda")
  expect_error(SVC_selection_control(n.lambda = 1.5), "n.lambda")
  expect_error(SVC_selection_control(method = "MBO", n.init = 2), "n.init")
  expect_error(SVC_selection_control(n.iter = NA_real_), "n.iter")
  expect_error(SVC_selection_control(CD.conv = list(N = 2)), "CD.conv")
  expect_error(
    SVC_selection_control(CD.conv = list(N = 2, delta = 0, logLik = TRUE)),
    "CD.conv\\$delta"
  )
  expect_error(
    SVC_selection_control(CD.conv = list(N = 2, delta = 1e-6, logLik = NA)),
    "CD.conv\\$logLik"
  )
  expect_error(SVC_selection_control(hessian = c(TRUE, FALSE)), "hessian")
  expect_error(SVC_selection_control(adaptive = NA), "adaptive")
  expect_error(SVC_selection_control(parallel = TRUE), "parallel")
  expect_error(SVC_selection_control(optim.args = TRUE), "optim.args")
})

test_that("SVC_selection validates its objective and covariance parameters", {
  expect_error(SVC_selection(1, 1), "obj.fun.*SVC_obj_fun")
  expect_error(
    SVC_selection(structure(list(), class = "SVC_obj_fun"), 1),
    "malformed"
  )

  obj.fun <- selection_objective()
  expect_error(
    SVC_selection(obj.fun, mle.par = c(1, 1)),
    "mle.par.*length 3"
  )
  expect_error(
    SVC_selection(obj.fun, mle.par = c(1, 1, NA_real_)),
    "mle.par.*finite"
  )
  expect_error(SVC_selection(obj.fun, mle.par = c(1, 1, 1), control = TRUE),
               "control")
})

test_that("SVC_selection accepts a valid extracted objective", {
  obj.fun <- selection_objective()
  fit <- SVC_selection(
    obj.fun,
    mle.par = c(0.5, 0.5, 0.5),
    control = SVC_selection_control(
      r.lambda = c(1e-4, 1e-3),
      n.lambda = 1L,
      CD.conv = list(N = 1L, delta = 1e-4, logLik = TRUE),
      optim.args = list(
        lower = rep(1e-6, 3L),
        upper = rep(Inf, 3L),
        control = list(maxit = 5L)
      )
    )
  )

  expect_s3_class(fit, "SVC_selection")
  expect_length(fit$lambda, 2L)
  expect_true(is.finite(fit$PMLE_pars$final.IC))
})

test_that("SVC_selection accepts an objective extracted from a formula", {
  skip_if_not_installed("glmnet")
  data(SVCdata)
  id <- seq_len(12)
  df <- data.frame(y = SVCdata$y[id], x = SVCdata$X[id, 2])
  obj <- SVC_mle(
    y ~ x, data = df, RE_formula = ~ 1, locs = SVCdata$locs[id],
    control = SVC_mle_control(extract_fun = TRUE, profileLik = TRUE)
  )
  fit <- SVC_selection(
    obj,
    mle.par = c(0.5, 0.5, 0.5),
    control = SVC_selection_control(
      r.lambda = c(1e-4, 1e-3),
      n.lambda = 1L,
      CD.conv = list(N = 1L, delta = 1e-4, logLik = TRUE),
      optim.args = list(
        lower = rep(1e-6, 3L),
        upper = rep(Inf, 3L),
        control = list(maxit = 5L)
      )
    )
  )

  expect_s3_class(fit, "SVC_selection")
  expect_true(is.finite(fit$PMLE_pars$final.IC))
})
