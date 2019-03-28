# holds objective functions to be optimized, i.e. negative log-likelihood of SVC-Models


nLL <- function(x, cov_func, outer.X, y, X, X.fixed) {

  p <- length(outer.X)
  n <- length(y)

  # compute covariance matrices
  Sigma <- Sigma_y(x, p, cov_func, outer.X)


  # calculate Cholesky-Decompisition
  cholS <- spam::chol.spam(Sigma)

  # get beta
  if (is.null(X.fixed)) {
    beta <- x[1 + 2*p + 1:p]

    resid <- y - X %*% beta
  } else {
    p.fixed <- ncol(X.fixed)
    beta <- x[1 + 2*p + 1:(p + p.fixed)]

    resid <- y - cbind(X, X.fixed) %*% beta
  }

  return(n * log(2 * pi) +
           2 * c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus) +
           as.numeric(crossprod(spam::forwardsolve(cholS, resid,
                                                   transpose = TRUE,
                                                   upper.tri = TRUE))))

}
