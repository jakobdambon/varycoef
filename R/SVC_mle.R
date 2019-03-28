## -----------------------------------------------------------------------------
## In this script, one finds every function directly related to estimating
## and predicting SVC using our proposed MLE.
## -----------------------------------------------------------------------------



## ---- help function to give back correct covariance function ----
MLE.cov.func <- function(cov.name) {
  cov.func = switch(cov.name,
                    "exp" = spam::cov.exp,
                    # "sph" = spam::cov.sph,
                    stop("SVC.cov argument not defined."))
}




#' @title MLE of SVC model
#'
#' @description Calls MLE of the SVC model defined as:
#'
#' \deqn{y(s) = x^{(1)}(s)\beta_1(s) + ... + x^{(p)}(s)\beta_p(s) + x^{(p+1)}\beta_{p+1} + ... + x^{(p+p.fix)}\beta_{p+p.fix} + \epsilon(s)}
#'
#' where:
#' \itemize{
#'   \item p is the number of SVC
#'   \item \eqn{\beta_j} are distinctly defined Gaussian Random Fields, j = 1, ..., p
#'   \item p.fix is the number of fixed coefficients, which is optional
#'   \item \eqn{\epsilon} is the nugget effect
#' }
#'
#' The MLE is done by calling the function optim
#' @param y         numeric response vector of dimension n.
#' @param X         matrix of covariates of dimension n x p. Intercept has to be added manually.
#' @param locs      matrix of locations of dimension n X 2. May contain multiple observations at single location which (may) cause a permutation of \code{y}, \code{X}, \code{X.fixed} and \code{locs}.
#' @param init      numeric. Initial values for optimization procedure. The vector consists of p-times (alternating) scale and variance, the nugget variance and the p + p.fix mean effects
#' @param control   list of control paramaters, see Details
#' @param X.fixed   Do not use this argument. Optional matrix of covariates with fixed effects, i.e. non-SVC, of dimension n x p.fixed
#' @param ns        Do not use this argument.
#' @param ...       futher arguments for method
#'
#'
#' @return Object of class \code{SVC_mle}
#'
#' @details The \code{control} list has the following components
#' \itemize{
#'   \item \code{cov.name}, name of the covariance function defining the covariance matrix of the GRF. Currently, only \code{"exp"} for the exponential is supported.
#'   \item \code{tapering}, if \code{NULL}, no tapering is applied. If a scalar is given, covariance tapering with this taper range is applied.
#'   \item \code{cl}, cluster for parallelization. Currently not supported.
#'   \item \code{X.scale}, if \code{TRUE}, the covariates are being standardized. Currently not supported.
#' }
#'
#' @import spam
#' @import methods
#' @export
SVC_mle <- function(y, X, locs, init,
                    control = list(cov.name = "exp",
                                   tapering = NULL,
                                   cl = NULL,
                                   X.scale = FALSE),
                    X.fixed = NULL, ns = NULL, ...) {

  # check for multiple observations at locations
  if (nrow(unique(locs)) < nrow(locs)) {
    warning("Multiple Observations at single location detected.\nPermuting Observations!")
    # ordering by location
    u.locs <- unique(locs)
    ch.locs <- apply(locs, 1, paste0, collapse = "x")
    u.ch.locs <- unique(sort(ch.locs))
    ord.by.locs <- order(ch.locs)

    ns <- as.numeric(table(ch.locs)[u.ch.locs])

    J <- spam::diag.spam(nrow(locs))
    J@colindices <- unlist(mapply(rep, times = ns, x = 1:nrow(u.locs)))
    J@dimension[2] <- nrow(u.locs)

    X.tilde <- sapply(1:ncol(X), function(j){
      spam::diag.spam(spam::crossprod.spam(J, spam::diag.spam(X[ord.by.locs, j]))%*%J)
    })

    X.fixed.tilde <- if (is.null(X.fixed)) {
      NULL
    } else {
      sapply(1:ncol(X.fixed), function(j){
        spam::diag.spam(spam::crossprod.spam(J, spam::diag.spam(X.fixed[ord.by.locs, j]))%*%J)
      })
    }

    return(SVC_mle(y = spam::crossprod.spam(J, y[ord.by.locs]),
                   X = X.tilde,
                   locs = u.locs[order(u.ch.locs)],
                   init = init,
                   control = control,
                   cl = cl,
                   X.fixed = X.fixed.tilde,
                   ns = ns, ...))
  }


  # scaling?
  if (control$X.scale) {
    scale.pars <- apply(X, 2, function(x) {c(m = mean(x), s = sd(x))})
    # scaled X = Z
    Z <- apply(1:ncol(X), 2, function(j) {
      m.s <- scale.pars[, j]
      if (m.s[2] != 0) {
        return((X[, j]-m.s[1])/m.s[2])
      } else {
        return(X[, j])
      }
    })

    # scale X.fixed if necessary
    if (!is.null(X.fixed)) {
      scale.pars.fixed <- apply(X.fixed, 2, function(x) {c(m = mean(x), s = sd(x))})
      # scaled X = Z
      Z.fixed <- apply(1:ncol(X.fixed), 2, function(j) {
        m.s <- scale.pars.fixed[, j]
        if (m.s[2] != 0) {
          return((X.fixed[, j]-m.s[1])/m.s[2])
        } else {
          return(X.fixed[, j])
        }
      })
    } else {
      Z.fixed <- scale.pars.fixed <- NULL
    }

    control$scale.pars <- list(scale.pars = scale.pars,
                               scale.pars.fixed = scale.pars.fixed)


    # call with scaled matrices
    return(SVC_mle(y = y,
                   X = Z,
                   locs = locs,
                   init = init,
                   control = control,
                   cl = cl,
                   X.fixed = Z.fixed,
                   ns = ns,
                   ...))

  } # end scaling


  # define distance matrix
  if (is.null(control$tapering)) {
    d <- spam::as.spam(dist(locs))
  } else {
    d <- spam::nearest.dist(locs, delta = control$tapering)
  }


  # get covariance function
  raw.cov.func <- MLE.cov.func(control$cov.name)

  cov.func <- list(
    # covariance function
    cov.func =
      if(is.null(control$taper)) {
        # without tapering
        function(x) raw.cov.func(d, x)
      } else {
        # with tapering
        function(x)
          raw.cov.func(d, x) *
          spam::cov.wend1(d, c(control$taper.r, 1, 0))
      },
    # number of observations at single location (needed for nugget)
    ns = ns)

  # get outer matrices of observations
  outer.X <- lapply(1:ncol(X), function(j) X[, j]%o%X[, j])

  # call optimization
  if (is.null(control$cl)) {

  } else {
    stop("Parallelization not yet implemented.")
  }


  optim.output <- optim(par = init,
                        fn = nLL,
                          # arguments of nLL
                          cov_func = cov.func,
                          outer.X  = outer.X,
                          y        = y,
                          X        = X,
                          X.fixed  = X.fixed,
                        method = "L-BFGS-B",
                        ...)

  result <- list(optim.output = optim.output,
                 call.args = list(y = y,
                                  X = X,
                                  locs = locs,
                                  init = init,
                                  control = control,
                                  X.fixed = X.fixed,
                                  ns = ns, ...),
                 comp.args = list(outer.X = outer.X))

  class(result) <- "SVC_mle"
  return(result)
}




#' Prediction of SVC (and response variable)
#'
#' @param object   output of \code{SVC_mle}
#' @param newlocs matrix of dimension n' x 2. These are the new locations the SVCs are predicted for. If \code{NULL}, the locations from the \code{SVC_mle} (i.e. \code{locs}) are considered.
#' @param newX    optional matrix of dimension n' x p. If provided, besides the predicted SVC, the function also returns the predicted response variable.
#' @param ...     further arguments
#' @return returns a data frame of n' rows and with columns
#' \itemize{
#'   \item \code{SVC_1, ..., SVC_p}, i.e. the predicted SVC at locations \code{newlocs}
#'   \item \code{y.pred}, if \code{newX} is provided
#'   \item \code{loc_x, loc_y}, the locations of the predictions
#' }
#' @import spam
#' @importFrom fields rdist
#' @export
predict.SVC_mle <- function(object, newlocs = NULL, newX = NULL, ...) {

  hyper.par <- object$optim.output$par

  p <- length(object$comp.args$outer.X)
  n <- length(object$call.args$y)

  # if no new locations are given,
  # predict for training data
  if (is.null(newlocs)) {
    # compute untapered distance matrix
    newlocs <- object$call.args$locs
    d <- dd <- as.matrix(dist(newlocs))
    d[base::upper.tri(dd, diag = TRUE)] <- 0
    n.new <- n
  } else {
    d <- as.matrix(dist(object$call.args$locs))
    d[base::upper.tri(d)] <- 0
    dd <- fields::rdist(newlocs, object$call.args$locs)
    n.new <- nrow(newlocs)
  }

  # covariance function (not tapered)
  raw.cf <- MLE.cov.func(object$call.args$control$cov.name)

  # covariance y
  cf <- function(x) raw.cf(spam::as.spam(d), x)
  # cross-covariance (newlocs and locs)
  cf_dd <- function(x) raw.cf(dd, x)

  cov_y <- Sigma_y(x = hyper.par,
                   p = p,
                   cov_func = list(cov.func = cf, ns = object$call.args$ns),
                   outer.X = object$comp.args$outer.X)


  # cross-covariance beta' y
  cov_b_y <- Sigma_b_y(x = hyper.par,
                       cov.func = cf_dd,
                       X = as.matrix(object$call.args$X),
                       n.new = n.new)



  eff <- rep(hyper.par[2*p + 1 + 1:p], each = n.new) +
    cov_b_y %*% spam::solve.spam(cov_y) %*%
    (object$call.args$y - object$call.args$X %*% hyper.par[2*p + 1 + 1:p])

  eff <- matrix(eff, ncol = length(object$comp.args$outer.X))

  if (!is.null(newX)) {

    y.pred <- apply(as.matrix(newX) * eff, 1, sum)

    out <- as.data.frame(cbind(eff, y.pred, newlocs))
    colnames(out) <- c(paste0("SVC_", 1:ncol(eff)), "y.pred", "loc_x", "loc_y")
    return(out)
  } else {
    out <- as.data.frame(cbind(eff, newlocs))
    colnames(out) <- c(paste0("SVC_", 1:ncol(eff)), "loc_x", "loc_y")
    return(out)
  }


}


