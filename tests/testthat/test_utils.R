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