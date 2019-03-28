Sigma_y <- function(x, p, cov_func, outer.X) {
  n <- nrow(outer.X[[1]])

  Sigma <- lapply(1:(p+1), function(j) {
    if (j == (p+1)) {
      hyper.par <- x[2*p+1]
      if (is.null(cov_func$ns)){
        return(spam::diag.spam(rep(hyper.par, n)))
      } else {
        return(spam::diag.spam(hyper.par*cov_func$ns))
      }
    } else {
      hyper.par <- c(x[2*(j-1) + 1:2], 0)
      Cov <- do.call(cov_func$cov.func, list(hyper.par))

      # is needed since as.spam.dist does not return 0 on diagonal
      spam::diag.spam(Cov) <- x[2*j]

      return(Cov * outer.X[[j]])
    }

  })

  options(spam.trivalues = TRUE)

  # add lower tri. cov-matrices up and mirror them to get full cov-matrix
  Sigma <- Reduce('+', Sigma)
  spam::lower.tri.spam(Sigma) + spam::t.spam(Sigma)
}


Sigma_b_y <- function(x, cov.func, X, n.new) {
  n <- nrow(X)
  
  cov <- lapply(1:ncol(X),
                function(j) {
                  # cross-covariances of Sigma_b_y
                  cov.func(c(x[2*(j-1) + 1:2], 0)) * 
                    matrix(rep(X[, j], each = n.new), ncol = n)
                })
  # binding to one matrix
  Reduce(rbind, cov)
}


Sigma_y_y <- function(x, cov.func, X, newX) {
  n <- nrow(X)
  
  cov <- lapply(1:ncol(X),
                function(j) {
                  # cross-covariances of Sigma_b_y
                  cov.func(c(x[2*(j-1) + 1:2], 0)) * 
                    (newX[, j] %o% X[, j])
                })
  # binding to one matrix
  Reduce('+', cov)
}
