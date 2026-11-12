tr <- function(A) {
  # computes the trace of a matrix
  sum(diag(A))
}


eff_dof <- function(cov.par, X, cov_func, outer.W, taper) {

  n <- nrow(X)
  p <- length(outer.W)
  nugget_var <- cov.par[length(cov.par)]
  covariance <- Sigma_y(cov.par, cov_func, outer.W, taper)

  inverse_covariance <- solve(covariance)

  design_t_inverse <- crossprod(X, inverse_covariance)

  # trace of hat matrix (equation 15 in Mueller et al. 2013, Stat. Sci.)
  as.numeric(nugget_var * tr(solve(design_t_inverse %*% X) %*%
                               design_t_inverse %*% inverse_covariance %*% X) +
               n - nugget_var * tr(inverse_covariance))
}

tr_Sigma <- function(cov.par, X, cov_func, outer.W, taper) {

  n <- nrow(X)
  nugget_var <- tail(cov.par, 1)
  p <- length(outer.W)


  tr(Sigma_y(cov.par, cov_func, outer.W, taper)) - nugget_var*n
}

