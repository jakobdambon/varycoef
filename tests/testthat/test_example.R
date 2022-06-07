test_that("structure and dimensions of sample_SVCdata are correct", {
  n <- 10L
  
  set.seed(123)
  # SVC parameters
  df.pars <- data.frame(
    var = c(2, 0, 1),
    scale = c(3, 1, 1),
    mean = c(1, 2, 0)
  )
  # nugget standard deviation
  tau <- 0.5
  
  # sample locations
  s <- sort(runif(n, min = 0, max = 10))
  SVCdata <- sample_SVCdata(
    df.pars = df.pars, nugget.sd = tau, locs = s, cov.name = "mat32"
  )
  
  expect_type(SVCdata, "list")
  expect_identical(names(SVCdata), c("y", "X", "beta", "eps", "locs", "true_pars"))
  expect_length(SVCdata$y, n)
  expect_identical(dim(SVCdata$X), c(n, 3L))
  expect_identical(dim(SVCdata$beta), c(n, 3L))
  expect_length(SVCdata$eps, n)
  expect_type(SVCdata$true_pars, "list")
  expect_identical(dim(SVCdata$true_pars), c(4L, 3L))
  expect_identical(SVCdata$true_pars[1:3, ], df.pars)
  expect_identical(SVCdata$true_pars[4, 1], tau^2)
  expect_identical(SVCdata$true_pars[4, 2], NA_real_)
  expect_identical(SVCdata$true_pars[4, 3], NA_real_)
})

test_that("all possible four types of SVCs are possible", {
  n <- 10L
  
  set.seed(123)
  # SVC parameters
  df.pars <- data.frame(
    var = c(2, 1, 0, 0),
    scale = c(3, 1, 1, 1),
    mean = c(10, 0, 2, 0)
  )
  # nugget standard deviation
  tau <- 0.5
  
  # sample locations
  s <- sort(runif(n, min = 0, max = 10))
  SVCdata <- sample_SVCdata(
    df.pars = df.pars, nugget.sd = tau, locs = s, cov.name = "mat32"
  )
  
  # first SVC (FE + RE):
  # mean greater 0
  expect_true(all(SVCdata$beta[, 1] > 1))
  # variation
  expect_true(min(SVCdata$beta[, 1]) != max(SVCdata$beta[, 1]))
  
  # second SVC (only RE): 
  # variation
  expect_true(min(SVCdata$beta[, 2]) != max(SVCdata$beta[, 2]))
  
  # third SVC (only FE): 
  # mean equal to 2
  expect_true(all(SVCdata$beta[, 3] == 2))

  # fourth SVC (neither FE nor RE): 
  # mean equal to 2
  expect_true(all(SVCdata$beta[, 4] == 0))
})