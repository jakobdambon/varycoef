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
