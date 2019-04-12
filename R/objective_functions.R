# holds objective functions to be optimized, i.e. negative log-likelihood of SVC-Models


nLL <- function(x, cov_func, outer.W, y, X, W, taper = NULL) {


  pW <- ncol(W)
  pX <- ncol(X)
  n <- length(y)

  # compute covariance matrices
  Sigma <- Sigma_y(x, pW, cov_func, outer.W, taper = taper)


  # calculate Cholesky-Decompisition
  cholS <- spam::chol.spam(Sigma)

  # get beta
  beta <- x[1 + 2*pW + 1:pX]

  resid <- y - X %*% beta


  return(n * log(2 * pi) +
           2 * c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus) +
           as.numeric(crossprod(spam::forwardsolve(cholS, resid,
                                                   transpose = TRUE,
                                                   upper.tri = TRUE))))


}
