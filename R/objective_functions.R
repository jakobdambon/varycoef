# holds objective functions to be optimized, i.e. negative log-likelihood of SVC-Models

profile.n2LL <- function(x, cov_func, outer.W, y, X, W, mean.est,
                         taper = NULL, envir = NULL, pc.dens = NULL) {

  pW <- ncol(W)
  pX <- ncol(X)
  n <- length(y)

  # compute covariance matrices
  Sigma <- Sigma_y(x, pW, cov_func, outer.W, taper = taper)

  # calculate Cholesky-Decompisition
  cholS <- spam::chol.spam(Sigma)


  ## profile LL
  # compute mu(theta)

  # inverse of Sigma
  I.Sigma <- solve.spam(Sigma, Rstruct = cholS)

  # compute mu (using GLS) or take OLS mu, which is doe not change
  if (is.null(mean.est)) {
    # compute mu
    B <- crossprod.spam(X, I.Sigma)

    mu <- solve(B %*% X) %*% B %*% y
  } else {
    mu <- mean.est
  }


  # save past computations of mu
  if (!is.null(envir)) {
    envir$mu <- cbind(envir$mu, mu)
    envir$x  <- cbind(envir$x, matrix(x, ncol = 1))
  }

  res <- y - X %*% mu


  return(n * log(2 * pi) +
           2 * c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus) +
           as.numeric(crossprod(res, I.Sigma %*% res)) +
           pc_penalty(x, pW, pc.dens))

}


n2LL <- function(x, cov_func, outer.W, y, X, W,
                 taper = NULL, envir = NULL, pc.dens = NULL) {


  pW <- ncol(W)
  pX <- ncol(X)
  n <- length(y)

  # compute covariance matrices
  Sigma <- Sigma_y(x, pW, cov_func, outer.W, taper = taper)


  # calculate Cholesky-Decompisition
  cholS <- spam::chol.spam(Sigma)

  # get mu
  mu <- x[1 + 2*pW + 1:pX]

  res <- y - X %*% mu

  # save all x, mu
  if (!is.null(envir)) {
    envir$x  <- cbind(envir$x, matrix(x[1:(2*pW+1)], ncol = 1))
    envir$mu <- cbind(envir$mu, mu)
  }


  return(n * log(2 * pi) +
           2 * c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus) +
           as.numeric(crossprod(spam::forwardsolve(cholS, res,
                                                   transpose = TRUE,
                                                   upper.tri = TRUE)))+
           pc_penalty(x, pW, pc.dens))


}


pc_penalty <- function(x, pW, pc.dens) {

  if (is.null(pc.dens)) {
    return(0)
  } else {
    cov_vars <- x[1:(2*pW)]

    pc.priors <- sapply(1:pW, function(j) {
      pc.dens(cov_vars[2*(j-1) + 1:2])
    })

    return(sum(pc.priors))
  }
}

