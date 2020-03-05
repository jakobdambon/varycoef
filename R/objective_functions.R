# holds objective functions to be optimized, i.e. negative log-likelihood of SVC-Models

profile.n2LL <- function(x, cov_func, outer.W, y, X, W, mean.est,
                         taper = NULL, pc.dens = NULL) {

  pW <- ncol(W)
  pX <- ncol(X)
  n <- length(y)

  # compute covariance matrices
  Sigma <- spam::as.spam(Sigma_y(x, pW, cov_func, outer.W, taper = taper))

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

  res <- y - X %*% mu

  # n2LL as stated in paper Dambon et al. (2020) does not contain the
  # summand n*log(2 * pi), but it is needed to compute the actual LL
  return(n * log(2 * pi) +
           2 * c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus) +
           as.numeric(crossprod(res, I.Sigma %*% res)) +
           pc_penalty(x, pW, pc.dens))

}


n2LL <- function(x, cov_func, outer.W, y, X, W,
                 taper = NULL, pc.dens = NULL) {


  pW <- ncol(W)
  pX <- ncol(X)
  n <- length(y)

  # compute covariance matrices
  Sigma <- Sigma_y(x, pW, cov_func, outer.W, taper = taper)


  # calculate Cholesky-Decompisition
  cholS <- spam::chol.spam(spam::as.spam(Sigma))

  # get mu
  mu <- x[1 + 2*pW + 1:pX]

  res <- y - X %*% mu


  # n2LL as stated in paper Dambon et al. (2020) does not contain the
  # summand n*log(2 * pi), but it is needed to compute the actual LL
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

