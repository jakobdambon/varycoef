Sigma_y <- function(x, p, cov_func, outer.W, taper = NULL) {
  n <- nrow(outer.W[[1]])

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

      return(Cov * outer.W[[j]])
    }

  })

  options(spam.trivalues = TRUE)

  
  Sigma <- Reduce('+', Sigma)
  if (!is.null(taper)) {
    Sigma <- Sigma * taper
  }
  # add lower tri. cov-matrices up and mirror them to get full cov-matrix
  spam::lower.tri.spam(Sigma) + spam::t.spam(Sigma)
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


par.mapping <- function(p, ind.0.W, ind.0.X) {
  return(function(x) {
    out <- rep(NA, 3*p+1)
    jW <- 1
    jX <- 1
    for (j in 1:p) {
      if (j %in% ind.0.W) {
        out[2*(j-1)+1:2] <- c(1, 1)
      } else {
        out[2*(j-1)+1:2] <- x[2*(jW-1)+1:2]
        jW <- jW+1
      }

      if (j %in% ind.0.X) {
        out[2*p+1+j] <- 0
      } else {
        out[2*p+1+j] <- x[2*(p-length(ind.0.W))+1+jX]
        jX <- jX+1
      }
    }
    out[2*p+1] <- x[2*(p-length(ind.0.W))+1]
    return(out)
  })
}
