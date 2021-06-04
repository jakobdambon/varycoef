# holds objective functions to be optimized, i.e. negative log-likelihood of SVC-Models

n2LL <- function(
  x, cov_func, outer.W, y, X, W,
  mean.est = NULL,
  taper = NULL,
  pc.dens = NULL,
  Rstruct = NULL,
  profile = TRUE
) {

  q <- dim(W)[2]
  p <- dim(X)[2]
  n <- length(y)

  # compute covariance matrices
  Sigma <- Sigma_y(x[1:(2*q+1)], cov_func, outer.W, taper = taper)

  # calculate Cholesky-Decompisition
  # powerboost function
  # spam pivot check
  if (is.spam(Sigma)) {
    cholS <- spam::chol.spam(Sigma, Rstruct = Rstruct)
  } else {
    cholS <- chol(Sigma)
  }


  ## profile LL
  mu <- if (profile) {
    # compute mu(theta)...
    if (is.null(mean.est)) {
      # ...using GLS
      GLS_chol(cholS, X, y)
    } else {
      # ...or set by some constant
      mean.est
    }
  } else {
    # given directly by objective parameters
    x[1 + 2*q + 1:p]
  }

  res <- y - X %*% mu

  ## quadratic form of residual, i.e., t(res) %*% Sigma^-1 %*% res
  quad_res <- if (is.matrix(cholS)) {
    as.numeric(crossprod(solve(t(cholS), res)))
  } else {
    as.numeric(crossprod(spam::forwardsolve(cholS, res,
                                            transpose = TRUE,
                                            upper.tri = TRUE)))
  }

  # n2LL as stated in paper Dambon et al. (2020) does not contain the
  # summand n*log(2 * pi), but it is needed to compute the actual LL
  return(n * log(2 * pi) +
           2 * c(determinant(cholS)$modulus) +
           quad_res +
           pc_penalty(x, q, pc.dens))
}

#
# n2LL <- function(x, cov_func, outer.W, y, X, W,
#                  taper = NULL, pc.dens = NULL, Rstruct = NULL) {
#
#
#   pW <- ncol(W)
#   pX <- ncol(X)
#   n <- length(y)
#
#   # compute covariance matrices
#   Sigma <- Sigma_y(x, cov_func, outer.W, taper = taper)
#
#
#   if (is.spam(Sigma)) {
#     cholS <- spam::chol.spam(Sigma, Rstruct = Rstruct)
#   } else {
#     cholS <- chol(Sigma)
#   }
#
#   # calculate Cholesky-Decompisition
#   cholS <- spam::chol.spam(Sigma, Rstruct = Rstruct)
#
#   # get mu
#   mu <- x[1 + 2*pW + 1:pX]
#
#   res <- y - X %*% mu
#
#
#   quad_res <- if (is.matrix(cholS)) {
#     as.numeric(crossprod(solve(t(cholS), res)))
#   } else {
#     as.numeric(crossprod(spam::forwardsolve(cholS, res,
#                                             transpose = TRUE,
#                                             upper.tri = TRUE)))
#   }
#
#
#   # n2LL as stated in paper Dambon et al. (2020) does not contain the
#   # summand n*log(2 * pi), but it is needed to compute the actual LL
#   return(n * log(2 * pi) +
#            2 * c(determinant(cholS)$modulus) +
#            quad_res +
#            pc_penalty(x, pW, pc.dens))
# }


pc_penalty <- function(x, q, pc.dens) {

  if (is.null(pc.dens)) {
    return(0)
  } else {
    cov_vars <- x[1:(2*q)]

    pc.priors <- sapply(1:q, function(j) {
      pc.dens(cov_vars[2*(j-1) + 1:2])
    })

    return(sum(pc.priors))
  }
}
