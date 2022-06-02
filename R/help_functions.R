

#' GLS Estimate using Cholesky Factor
#'
#' Computes the GLS estimate using the formula:
#' \deqn{\mu_{GLS} = (X^\top \Sigma^{-1} X)^{-1}X^\top \Sigma^{-1} y.}
#' The computation is done depending on the input class of the Cholesky factor
#' \code{R}. It relies on the classical \code{\link[base]{solve}} or on
#' using \code{forwardsolve} and \code{backsolve} functions of package
#' \code{spam}, see \code{\link[spam]{solve}}. This is much faster than
#' computing the inverse of \eqn{\Sigma}, especially since we have to compute
#' the Cholesky decomposition of \eqn{\Sigma} either way.
#'
#' @param R (\code{spam.chol.NgPeyton} or \code{matrix(n, n)}) \cr Cholesky factor of
#' the covariance matrix \eqn{\Sigma}. If covariance tapering and sparse
#' matrices are used, then the input is of class \code{spam.chol.NgPeyton}.
#' Otherwise, \code{R} is the output of a standard \code{\link[base]{chol}},
#' i.e., a simple \code{matrix}
#' @param X (\code{matrix(n, p)}) \cr Data / design matrix.
#' @param y (\code{numeric(n)}) \cr Response vector
#'
#' @return A \code{numeric(p)} vector, i.e., the mean effects.
#' @author Jakob Dambon
#'
#' @export
#'
#' @examples
#' # generate data
#' n <- 10
#' X <- cbind(1, 20+1:n)
#' y <- rnorm(n)
#' A <- matrix(runif(n^2)*2-1, ncol=n)
#' Sigma <- t(A) %*% A
#' # two possibilities
#' ## using standard Cholesky decomposition
#' R_mat <- chol(Sigma); str(R_mat)
#' mu_mat <- GLS_chol(R_mat, X, y)
#' ## using spam
#' R_spam <- chol(as.spam(Sigma)); str(R_spam)
#' mu_spam <- GLS_chol(R_spam, X, y)
#' # should be identical to the following
#' mu <- solve(crossprod(X, solve(Sigma, X))) %*%
#'       crossprod(X, solve(Sigma, y))
#' ## check
#' abs(mu - mu_mat)
#' abs(mu - mu_spam)
GLS_chol <- function(R, X, y) UseMethod("GLS_chol")

#' @rdname GLS_chol
#' @export
GLS_chol.spam.chol.NgPeyton <- function(R, X, y) {
  # (X^T * Sigma^-1 * X)^-1
  solve(
    crossprod(spam::forwardsolve(R, X))
  ) %*%
    # (X^T * Sigma^-1 * y)
    crossprod(X, spam::backsolve(R, spam::forwardsolve(R, y)))
}

#' @rdname GLS_chol
#' @export
GLS_chol.matrix <- function(R, X, y) {
  RiX <- solve(t(R), X)
  # (X^T * Sigma^-1 * X)^-1
  solve(
    crossprod(RiX)
  ) %*%
    # (X^T * Sigma^-1 * y)
    crossprod(RiX, solve(t(R), y))
}


#' Covariance Matrix of GP-based SVC Model
#'
#' Builds the covariance matrix of \eqn{y} (p. 6, Dambon et al. (2021)
#'\doi{10.1016/j.spasta.2020.100470}) for a given set of covariance
#' parameters and other, pre-defined objects (like the outer-products,
#' covariance function, and, possibly, a taper matrix).
#'
#' @param x         (\code{numeric(2q+1)}) \cr Non negative vector containing
#' the covariance parameters in the following order: \eqn{\rho_1, \sigma_1^2,
#' ..., \rho_q, \sigma_q^2 , \tau^2}. Note that the odd entries, i.e., the
#' ranges and the nugget variance, have to be greater than 0, otherwise the
#' covariance matrix is not well-defined (singularities or not-invertible).
#' @param cov_func  (\code{function}) \cr A covariance function that works on
#' the pre-defined distance matrix \code{d}. It takes a numeric vector as an
#' input, the first entry being the range, the second being the variance
#' (also called partial sill). Usually, it is defined as, e.g.:
#' \code{function(pars) spam::cov.exp(d, pars)} or any other covariance function
#' defined for two parameters.
#' @param outer.W   (\code{list(q)}) \cr A list of length \code{q} containing
#' the outer products of the random effect covariates in a lower triangular,
#' (possibly sparse) matrix. If tapering is applied, the list entries, i.e.,
#' the outer products have to be given as \code{\link[spam]{spam}} objects.
#' @param taper     (\code{NULL} or \code{spam}) \cr If covariance tapering is
#' applied, this argument contains the taper matrix, which is a
#' \code{\link[spam]{spam}} object. Otherwise, it is \code{NULL}.
#'
#' @return Returns a positive-definite covariance matrix y, which is needed in
#' the MLE. Specifically, a Cholesky Decomposition is applied on the covariance
#' matrix.
#'
#'
#' @author Jakob Dambon
#' @references Dambon, J. A., Sigrist, F., Furrer, R. (2021)
#'    \emph{Maximum likelihood estimation of spatially varying coefficient
#'    models for large data with an application to real estate price prediction},
#'    Spatial Statistics \doi{10.1016/j.spasta.2020.100470}
#'
#'
#' @examples
#' # locations
#' locs <- 1:6
#' # random effects covariates
#' W <- cbind(rep(1, 6), 5:10)
#' # distance matrix with and without tapering
#' d <- as.matrix(dist(locs))
#' # distance matrix with and without tapering
#' tap_dist <- 2
#' d_tap <- spam::nearest.dist(locs, delta = tap_dist)
#' # call without tapering
#' (Sy <- varycoef:::Sigma_y(
#'   x = rep(0.5, 5),
#'   cov_func = function(x) spam::cov.exp(d, x),
#'   outer.W = lapply(1:ncol(W), function(k) W[, k] %o% W[, k])
#' ))
#' str(Sy)
#' # call with tapering
#' (Sy_tap <- varycoef:::Sigma_y(
#'   x = rep(0.5, 5),
#'   cov_func = function(x) spam::cov.exp(d_tap, x),
#'   outer.W = lapply(1:ncol(W), function(k)
#'     spam::as.spam((W[, k] %o% W[, k]) * (d_tap<=tap_dist))
#'   ),
#'   taper = spam::cov.wend1(d_tap, c(tap_dist, 1, 0))
#' ))
#' str(Sy_tap)
#' # difference between tapered and untapered covariance matrices
#' Sy-Sy_tap
Sigma_y <- function(x, cov_func, outer.W, taper = NULL) {
  n <- nrow(outer.W[[1]])
  q <- length(outer.W)

  if (is.null(taper)) {
    # with no tapering computations are done on matrix objects
    Sigma <- matrix(0, nrow = n, ncol = n)

    for (k in 1:q) {
      # first argument: range, second argument: variance / sill
      Cov <- cov_func(x[2*(k-1) + 1:2])

      Sigma <- Sigma + (
        Cov * outer.W[[k]]
      )
    }

    nug <- if (n == 1) {
      x[2*q+1]
    } else {
      diag(rep(x[2*q+1], n))
    }

    return(Sigma + nug)
  } else {
    # With tapering computations are done on spam objects.
    # Specifically, due to their fixed structure and since we are only
    # pair-wise adding and multiplying, on the spam entries themselves

    stopifnot(
      all(sapply(outer.W, is.spam))
    )

    # construct a sparse matrix with 0 values as future entries
    # for k = 1
    Sigma <- outer.W[[1]] * cov_func(x[1:2])

    # if q > 1, build covariance matrix using components of other GPs
    if (q > 1) {
      for (k in 2:q) {
        Cov <- do.call(cov_func, list(c(x[2*(k-1) + 1:2], 0)))

        Sigma <- Sigma + (Cov * outer.W[[k]])
      }
    }

    options(spam.trivalues = TRUE)

    nug <- if (n == 1) {
      x[2*q+1]
    } else {
      spam::diag.spam(rep(x[2*q+1], n))
    }

    # Sigma <- Sigma * taper
    # add lower tri. cov-matrices up and mirror them to get full cov-matrix
    # due to spam::nearest.dist design

    return(spam::lower.tri.spam(Sigma) +
             spam::t.spam(Sigma) +
             nug)
  }
}




Sigma_b_y <- function(x, cov.func, W, n.new) {
  n <- nrow(W)

  cov <- lapply(1:ncol(W),
                function(j) {
                  # cross-covariances of Sigma_b_y
                  cov.func(c(x[2*(j-1) + 1:2])) *
                    matrix(rep(W[, j], each = n.new), ncol = n)
                })
  # binding to one matrix
  Reduce(rbind, cov)
}


Sigma_y_y <- function(x, cov.func, X, newX) {

  p <- ncol(X)

  Sigma <- matrix(0, ncol = nrow(X), nrow = nrow(newX))

  for (j in 1:p) {
    Cov <- cov.func(c(x[2*(j-1) + 1:2], 0))

    Sigma <- Sigma + (
      Cov * (newX[, j] %o% X[, j])
    )
  }


  return(Sigma)

}


#' Check Lower Bound of Covariance Parameters
#'
#' Ensures that the covariance parameters define a positive definite covariance
#' matrix. It takes the vector
#' \eqn{(\rho_1, \sigma^2_1, ..., \rho_q, \sigma^2_q, \tau^2)} and checks if
#' all \eqn{\rho_k>0}, all \eqn{\sigma_k^2>=0}, and \eqn{\tau^2>0}.
#' @param cv (\code{numeric(2*q+1)}) \cr Covariance vector of SVC model.
#' @param q  (\code{numeric(1)}) \cr Integer indicating the number of SVCs.
#'
#' @return \code{logical(1)} with \code{TRUE} if all conditions above are
#' fulfilled.
#' @export
#'
#' @examples
#' # first one is true, all other are false
#' check_cov_lower(c(0.1, 0, 0.2,  1, 0.2), q = 2)
#' check_cov_lower(c(0  , 0, 0.2,  1, 0.2), q = 2)
#' check_cov_lower(c(0.1, 0, 0.2,  1, 0  ), q = 2)
#' check_cov_lower(c(0.1, 0, 0.2, -1, 0  ), q = 2)
check_cov_lower <- function(cv, q) {
  # check range and nugget variance parameters
  l_rp <- all(cv[1+2*(0:q)]>0)
  # check SVC variances
  l_vp <- all(cv[2*(1:q)] >= 0)
  return(l_rp & l_vp)
}


#' Setting of Optimization Bounds and Initial Values
#'
#' Sets bounds and initial values for \code{\link[stats]{optim}} by
#' extracting potentially given values from \code{\link{SVC_mle_control}} and
#' checking them, or calculating them from given data. See Details.
#'
#' @param control  (\code{\link{SVC_mle_control}} output, i.e. \code{list})
#' @param p        (\code{numeric(1)}) \cr Number of fixed effects
#' @param q        (\code{numeric(1)}) \cr Number of SVCs
#' @param id_obj   (\code{numeric(2*q+1+q)}) \cr Index vector to identify the
#' arguments of objective function.
#' @param med_dist (\code{numeric(1)}) \cr Median distance between observations
#' @param y_var    (\code{numeric(1)}) \cr Variance of response \code{y}
#' @param OLS_mu   (\code{numeric(p)}) \cr Coefficient estimates of ordinary
#' least squares (OLS).
#'
#' @details If values are not provided, then they are set in the following way.
#'    Let \eqn{d} be the median distance \code{med_dist}, let \eqn{s^2_y} be
#'    the variance of the response \code{y_var}, and let \eqn{b_j} be the OLS
#'    coefficients of the linear model. The computed values are given in the
#'    table below.
#'
#'    | Parameter    | Lower bound   | Initial Value     | Upper Bound   |
#'    | ------------ | -------------:| -----------------:| -------------:|
#'    | Range        |  \eqn{d/1000} |         \eqn{d/4} |    \eqn{10 d} |
#'    | Variance     |       \eqn{0} | \eqn{s^2_y/(q+1)} | \eqn{10s^2_y} |
#'    | Nugget       | \eqn{10^{-6}} | \eqn{s^2_y/(q+1)} | \eqn{10s^2_y} |
#'    | Mean \eqn{j} |   \code{-Inf} |         \eqn{b_j} |    \code{Inf} |
#' @md
#'
#' @author Jakob Dambon
#'
#' @export
#'
#' @return A \code{list} with three entries: \code{lower}, \code{init},
#' and \code{upper}.
init_bounds_optim <- function(control, p, q, id_obj, med_dist, y_var, OLS_mu) {

  # lower bound for optim
  if (is.null(control$lower)) {
    lower <- if (control$profileLik) {
      c(rep(c(med_dist/1000, 0), q), 1e-6)
    } else {
      c(rep(c(med_dist/1000, 0), q), 1e-6, rep(-Inf, p))
    }
  } else {
    lower <- control$lower
    if (length(lower) != length(id_obj)) {
      stop("Lower boundary vector has wrong length. Check SVC_mle_control.")
    }
    if (!check_cov_lower(lower[1:(2*q+1)], q)) {
      stop("Lower boundary vector is not greater (or equal) than 0. Call ?check_cov_lower.")
    }
  }

  # upper bound for optim
  if (is.null(control$upper)) {
    upper <- if (control$profileLik) {
      c(rep(c(10*med_dist, 10*y_var), q), 10*y_var)
    } else {
      c(rep(c(10*med_dist, 10*y_var), q), 10*y_var, rep(Inf, p))
    }
  } else {
    upper <- control$upper
    if (length(upper) != length(id_obj)) {
      stop("Upper boundary vector has wrong length. Check SVC_mle_control.")
    }
    if (!check_cov_lower(upper[1:(2*q+1)], q)) {
      stop("Upper boundary vector is not greater (or equal) than 0. Call ?check_cov_lower.")
    }
    if (any(lower>upper)) {
      stop("Upper boundary vector smaller than lower boundary.")
    }
  }

  # init
  if (is.null(control$init)) {
    init <- if (control$profileLik) {
      c(rep(c(med_dist/4, y_var/(q+1)), q), y_var/(q+1))
    } else {
      c(rep(c(med_dist/4, y_var/(q+1)), q), y_var/(q+1), OLS_mu)
    }
  } else {
    init <- control$init
    if (length(init) != length(id_obj)) {
      stop("Initial vector has wrong length. Check SVC_mle_control.")
    }
    if (!check_cov_lower(init[1:(2*q+1)], q)) {
      stop("Initial value vector is not greater (or equal) than 0. Call ?check_cov_lower.")
    }
    if (!(all(lower <= init) & all(init <= upper))) {
      stop("Initial values do not lie between lower and upper boundarys.")
    }
  }

  return(list(lower = lower, init = init, upper = upper))
}



#' Preparation of Parameter Output
#'
#' Prepares and computes the ML estimates and their respective standard errors.
#' @param output_par  (\code{numeric}) \cr Found optimal value of
#' \code{\link[stats]{optim}}.
#' @param Sigma_final (\code{spam} or \code{matrix(n, n)}) \cr Covariance matrix
#' Sigma of SVC under final covariance parameters.
#' @param Rstruct     (\code{NULL} or \code{spam.chol.NgPeyton}) \cr If
#' covariance tapering is used, the Cholesky factor has been calculated
#' previously and can be used to efficiently update the Cholesky factor of
#' \code{Sigma_final}, which is an \code{spam} object.
#' @param profileLik  (\code{logical(1)}) \cr Indicates if optimization has been
#' conducted over full or profile likelihood.
#' @param X (\code{matrix(n, p)}) Design matrix
#' @param y (\code{numeric(p)}) Response vector
#' @param H (\code{NULL} or \code{matrix}) Hessian of MLE
#' @param q (\code{numeric(1)}) Number of SVC
#'
#' @return A \code{list} with two \code{data.frame}. Each contains the estimated
#' parameters with their standard errors of the fixed and random effects,
#' respectively.
#'
#' @importFrom methods is
prep_par_output <- function(output_par, Sigma_final, Rstruct, profileLik,
                            X, y, H, q) {
  p <- dim(as.matrix(X))[2]

  # get mean effects depending on likelihood optimization
  if (profileLik) {
    # calculate Cholesky-Decomposition
    if (is.spam(Sigma_final)) {
      cholS <- chol(Sigma_final, Rstruct = Rstruct)
    } else {
      cholS <- chol(Sigma_final)
    }

    mu <- GLS_chol(cholS, X, y)
  } else {
    mu <- output_par[2*q+1 + 1:p]
  }

  # get standard errors of parameters
  if (is.null(H)) {
    warning("MLE without Hessian. Cannot return standard errors of covariance parameters.")
    se_all <- rep(NA, length(output_par))
  } else {
    # divide by 2 due to (-2)*LL
    se_all <- try({sqrt(diag(solve(H/2)))}, silent = TRUE)

    # if no convergence, standard errors cannot be extracted
    if (methods::is(se_all, "try-error")) {
      warning("Could not invert Hessian.")
      se_all <- rep(NA, length(output_par))
    }
  }
  # on profile?
  if (profileLik) {
    se_RE <- se_all
    # compute variance covariance matrix for fixed effects

    # using GLS properties
    Sigma_FE <- solve(crossprod(X, solve(Sigma_final, X)))

    se_FE <- sqrt(diag(Sigma_FE))
  } else {
    se_RE <- se_all[1:(2*q+1)]
    se_FE <- se_all[2*q+1 + 1:p]
  }

  return(list(
    RE = data.frame(est = output_par[1:(2*q+1)], SE = se_RE),
    FE = data.frame(est = mu, SE = se_FE)
  ))
}

#
# own_dist <- function(
#   locs, newlocs = NULL, taper = NULL, method_list = NULL
# ) {
#   if (is.null(taper)) {
#     d <- as.matrix(
#       do.call(dist,
#         c(list(x = locs, diag = TRUE, upper = TRUE), method_list)))
#   } else {
#     d <- do.call(spam::nearest.dist,
#                  c(list(x = locs,
#                         delta  = control$tapering),
#                    control$dist))
#   }
# }


#' Computes (Cross-) Distances
#'
#' @param x     (\code{matrix}) \cr Matrix containing locations
#' @param y     (\code{NULL} or \code{matrix}) \cr If \code{NULL}, computes the
#'     distances between \code{x}. Otherwise, computes cross-distances, i.e.,
#'     pair-wise distances between rows of \code{x} and \code{y}.
#' @param taper (\code{NULL} or \code{numeric(1)}) \cr If \code{NULL}, all
#'     distances are considered. Otherwise, only distances shorter than
#'     \code{taper} are used. Hence the output will be a sparse matrix of type
#'     \code{\link[spam]{spam}}.
#' @param ...   Further arguments for either \code{\link[stats]{dist}} or
#'     \code{\link[spam]{nearest.dist}}.
#'
#' @return A \code{matrix} or \code{spam} object.
#' @importFrom spam nearest.dist
#' @importFrom stats dist
own_dist <- function(x, y = NULL, taper = NULL, ...) {

  d <- if (is.null(taper)) {
    # without tapering
    if (is.null(y)) {
      # no cross distances
      as.matrix(do.call(
        dist,
        c(list(x = x, diag = TRUE, upper = TRUE), ...)
      ))
    } else {
      # cross distances
      as.matrix(do.call(
        spam::nearest.dist,
        c(list(x = x, y = y, delta = 1e99), ...)
      ))
    }
  } else {
    # with tapering
    if (is.null(y)) {
      # no cross distances
      do.call(
        spam::nearest.dist,
        c(list(x = x, delta = taper), ...)
      )
    } else {
      # cross distances
      do.call(
        spam::nearest.dist,
        c(list(x = x, y = y, delta = taper), ...)
      )
    }
  }
  # return output
  d
}

get_taper <- function(cov.name, d, tapering) {
  switch(
    cov.name,
    "exp" = spam::cov.wend1(d, c(tapering, 1, 0)),
    "mat32" = spam::cov.wend1(d, c(tapering, 1, 0)),
    "mat52" = spam::cov.wend2(d, c(tapering, 1, 0))
  )
}
