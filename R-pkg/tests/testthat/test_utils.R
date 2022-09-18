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
  expect_identical(f(x, theta = cp), spam::cov.mat32(x, theta = cp))

  f <- MLE.cov.func("mat52")
  expect_identical(f(x, theta = cp), spam::cov.mat52(x, theta = cp))

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


test_that("GLS_chol utility function", {
  set.seed(123)
  n <- 1000
  D_spam <- spam::nearest.dist(
    seq(0, 1, length.out = n),
    delta = 0.2
  )
  
  cholS <- chol(as.matrix(exp(-D) + diag(rep(0.1, n))))
  cholS_spam <- spam::chol.spam(as.spam(S))
  
  X <- matrix(1:(2*n), ncol = 2)
  y <- matrix(rnorm(n), ncol = 1)
  
  # Check that we get the same: (X^T * Sigma^-1 * X)^-1 * (X^T * Sigma^-1 * y)
  expect_equal(
    solve( t(X) %*% solve(S) %*% X ) %*% (t(X) %*% solve(S) %*% y ),
    GLS_chol(cholS, X, y)
  )
  
  expect_equal(
    solve( t(X) %*% solve(S) %*% X ) %*% (t(X) %*% solve(S) %*% y ),
    GLS_chol(cholS_spam, X, y)
  )
  
  # redundant check
  expect_equal(
    GLS_chol(cholS_spam, X, y),
    GLS_chol(cholS, X, y)
  )
})
