Sigma_y <- function(x, p, cov_func, outer.W, taper = NULL) {
  n <- nrow(outer.W[[1]])

  Sigma <- matrix(0, nrow = n, ncol = n)

  for (j in 1:p) {
    Cov <- do.call(cov_func, list(c(x[2*(j-1) + 1:2], 0)))

    Sigma <- Sigma + (
      Cov * outer.W[[j]]
    )
  }

  options(spam.trivalues = TRUE)

  nug <- if (n == 1) {
    x[2*p+1]
  } else {
    spam::diag.spam(rep(x[2*p+1], n))
  }

  if (!is.null(taper)) {
    Sigma.tap <- taper * Sigma
    # add lower tri. cov-matrices up and mirror them to get full cov-matrix
    # due to spam::nearest.dist design
    return(spam::lower.tri.spam(Sigma.tap) +
             spam::t.spam(Sigma.tap) +
             nug)
  } else {
    return(Sigma + nug)
  }


}




Sigma_b_y <- function(x, cov.func, W, n.new) {
  n <- nrow(W)

  cov <- lapply(1:ncol(W),
                function(j) {
                  # cross-covariances of Sigma_b_y
                  cov.func(c(x[2*(j-1) + 1:2], 0)) *
                    matrix(rep(W[, j], each = n.new), ncol = n)
                })
  # binding to one matrix
  Reduce(rbind, cov)
}


Sigma_y_y <- function(x, cov.func, X, newX) {

  p <- ncol(X)

  Sigma <- matrix(0, ncol = nrow(X), nrow = nrow(newX))

  for (j in 1:p) {
    Cov <- do.call(cov.func, list(c(x[2*(j-1) + 1:2], 0)))

    Sigma <- Sigma + (
      Cov * (newX[, j] %o% X[, j])
    )
  }


  return(Sigma)

}

