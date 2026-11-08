test_that("formula functions works", {
  x <- as.character("y~X1 + X2")
  f1 <- as.formula(x)
  f2 <- ~X1 + X2
  
  expect_false(is.formula(x))
  expect_true(is.formula(f1))
  
  # need the as.character transformation as Environment attributes of the 
  # outputs are different
  expect_identical(
    as.character(varycoef:::drop_response(f1)), 
    as.character(f2)
  )
  expect_identical(
    as.character(varycoef:::drop_response(f2)), 
    as.character(f2)
  )
})


test_that("MLE.cov.func gets the right covariance functions", {
  set.seed(123)
  x <- c(0, runif(10), 1)
  # covariance parameters
  cp <- c(2/7, 4)
  
  f <- MLE.cov.func("exp")
  expect_identical(f(x, theta = cp), spam::cov.exp(x, theta = cp))
  
  f <- MLE.cov.func("mat32")
  expect_identical(f(x, theta = cp), spam::cov.mat(x, theta = c(cp, 3/2)))
  # covariance parameter length must be 2
  expect_error(f(x, theta = 1))
  expect_error(f(x, theta = rep(1, 3)))
  
  f <- MLE.cov.func("mat52")
  expect_identical(f(x, theta = cp), spam::cov.mat(x, theta = c(cp, 5/2)))
  # covariance parameter length must be 2
  expect_error(f(x, theta = 1))
  expect_error(f(x, theta = rep(1, 3)))
  
  f <- MLE.cov.func("sph")
  expect_identical(f(x, theta = cp), spam::cov.sph(x, theta = cp))
  
  f <- MLE.cov.func("wend1")
  expect_identical(f(x, theta = cp), spam::cov.wend1(x, theta = cp))
  
  f <- MLE.cov.func("wend2")
  expect_identical(f(x, theta = cp), spam::cov.wend2(x, theta = cp))
  
  # for function argument should return function
  f2 <- function(x) x^2
  f <- MLE.cov.func(f2)
  expect_identical(f(x), f2(x))
  
  expect_error(MLE.cov.func(2))
})

test_that("GLS_chol matrix solve matches the covariance formula", {
  set.seed(321)
  n <- 12L
  A <- matrix(rnorm(n * n), nrow = n)
  Sigma <- crossprod(A) + diag(n)
  R <- chol(Sigma)
  X <- cbind(1, rnorm(n), rnorm(n))
  y <- rnorm(n)

  expected <- solve(crossprod(X, solve(Sigma, X))) %*%
    crossprod(X, solve(Sigma, y))

  expect_equal(GLS_chol(R, X, y), expected, tolerance = 1e-10)
})

test_that("sparse distance initialization uses pairwise distances", {
  locs <- matrix(1:6, ncol = 1)
  distances <- spam::nearest.dist(locs, delta = 3)

  expect_equal(varycoef:::median_sparse_distance(distances), 2)
  expect_equal(
    varycoef:::median_sparse_distance(
      spam::nearest.dist(10 * locs, delta = 30)
    ),
    20
  )
})

test_that("sparse distance initialization requires a retained pair", {
  distances <- spam::nearest.dist(matrix(c(0, 10), ncol = 1), delta = 1)

  expect_error(
    varycoef:::median_sparse_distance(distances),
    "no non-zero pairwise distances are within the tapering range"
  )
})
