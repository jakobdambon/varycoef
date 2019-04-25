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



#' @title Set Parameters for \code{SVC_mle}
#'
#' @description Function to set up control parameters for \code{\link{SVC_mle}}
#'
#' @param cov.name name of the covariance function defining the covariance matrix of the GRF. Currently, only \code{"exp"} for the exponential is supported.
#' @param tapering if \code{NULL}, no tapering is applied. If a scalar is given, covariance tapering with this taper range is applied, for all GRF modelling the SVC.
#' @param cl       cluster for parallelization. Currently not supported.
#' @param Z.covars if \code{TRUE}, the covariates are being standardized. Not implemented.
#' @param init     numeric. Initial values for optimization procedure. The vector consists of p-times (alternating) scale and variance, the nugget variance and the p + p.fix mean effects
#' @param ...      further parameters yet to be implemented
#'
#' @return A list with which \code{\link{SVC_mle}} can be controlled
#' @export
#'
#' @seealso \code{\link{SVC_mle}}
#'
#' @examples
#' control <- SVC_mle.control()
#' control$init <- rep(0.3, 10)
SVC_mle.control <- function(cov.name = c("exp"),
                            tapering = NULL,
                            cl = NULL,
                            Z.covars = FALSE,
                            init = NULL, ...) {
  stopifnot(is.null(tapering) | (tapering>=0),
            is.logical(Z.covars))
  list(cov.name = match.arg(cov.name),
       tapering = tapering,
       cl = cl,
       Z.covars = Z.covars,
       init = init, ...)
}




#' @title MLE of SVC model
#'
#' @description Calls MLE of the SVC model defined as:
#'
#' \deqn{y(s) = X \mu + W \tilde \beta (s) + \epsilon(s)}
#'
#' where:
#' \itemize{
#'   \item y is the response (vector of length n)
#'   \item X is the data matrix for the fixed effects covariates
#'   \item \eqn{\mu} is the vetor containing the fixed effects
#'   \item W is the data matrix for the SVCs represented by zero mean GRF
#'   \item \eqn{\tilde \beta} are the SVCs represented by zero mean GRF
#'   \item \eqn{\epsilon} is the nugget effect
#' }
#'
#' The MLE is done by calling the function optim
#' @param y              numeric response vector of dimension n.
#' @param X              matrix of covariates of dimension n x pX. Intercept has to be added manually.
#' @param locs           matrix of locations of dimension n X 2. May contain multiple observations at single location which (may) cause a permutation of \code{y}, \code{X}, \code{X.fixed} and \code{locs}.
#' @param W              Optional matrix of covariates with fixed effects, i.e. non-SVC, of dimension n x pW
#' @param control        list of control paramaters, usually given by \code{\link{SVC_mle.control}}
#' @param ns             Do not use this argument.
#' @param optim.control  list of control arguments for optimization function, see Details in \code{\link{optim}}
#'
#'
#' @return Object of class \code{SVC_mle}
#'
#' @seealso \code{\link{SVC_mle}}, \code{\link{predict.SVC_mle}}
#'
#' @import spam
#' @import methods
#' @importFrom stats dist optim
#' @export
SVC_mle <- function(y, X, locs,
                    W = NULL,
                    control = SVC_mle.control(),
                    ns = NULL, 
                    optim.control = list()) {



  # check for arguments
  if (is.null(W)) {
    return(SVC_mle(y = y,
                   X = X,
                   locs = locs,
                   control = control,
                   W = X,
                   ns = NULL, 
                   optim.control = optim.control))
  }



  pW <- ncol(W)
  pX <- ncol(X)

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

    X.tilde <- sapply(1:pX, function(j){
      spam::diag.spam(spam::crossprod.spam(J, spam::diag.spam(X[ord.by.locs, j]))%*%J)
    })

    W.tilde <- if (is.null(W)) {
      NULL
    } else {
      sapply(1:pW, function(j){
        spam::diag.spam(spam::crossprod.spam(J, spam::diag.spam(W[ord.by.locs, j]))%*%J)
      })
    }

    return(SVC_mle(y = spam::crossprod.spam(J, y[ord.by.locs]),
                   X = X.tilde,
                   locs = u.locs[order(u.ch.locs), ],
                   control = control,
                   W = W.tilde,
                   ns = ns, 
                   optim.control = optim.control))
  }




  # standardizing?
  if (FALSE) {
    scale.pars.X <- apply(X, 2, function(x) {c(m = mean(x), s = sd(x))})
    # standardized X = Z
    Z <- sapply(1:pX, function(j) {
      m.s <- scale.pars.X[, j]
      if (m.s[2] != 0) {
        return((X[, j]-m.s[1])/m.s[2])
      } else {
        return(X[, j])
      }
    })


    scale.pars.W <- apply(W, 2, function(x) {c(m = mean(x), s = sd(x))})
    # standardized W = ZW
    ZW <- sapply(1:pW, function(j) {
      m.s <- scale.pars.W[, j]
      if (m.s[2] != 0) {
        return((W[, j]-m.s[1])/m.s[2])
      } else {
        return(W[, j])
      }
    })


    control$scale.pars <- list(scale.pars.X = scale.pars.X,
                               scale.pars.W = scale.pars.W)
    # overwrite standardization to avoid endless loop
    control$Z.covars <- FALSE
    # call with scaled matrices
    return(SVC_mle(y = y,
                   X = Z,
                   locs = locs,
                   control = control,
                   W = ZW,
                   ns = ns,
                   optim.control = optim.control))

  } # end scaling


  # define distance matrix
  if (is.null(control$tapering)) {
    d <- spam::as.spam(stats::dist(locs))
  } else {
    d <- spam::nearest.dist(locs, delta = control$tapering)
  }


  # get covariance function
  raw.cov.func <- MLE.cov.func(control$cov.name)

  cov.func <- list(
    # covariance function
    cov.func = function(x) raw.cov.func(d, x),
    # number of observations at single location (needed for nugget)
    ns = ns)

  # get outer matrices of observations
  outer.W <- lapply(1:pW, function(j) W[, j]%o%W[, j])


  # init
  if (is.null(control$init)) {
    init <- c(rep(0.3, 2*pW+1),
              rep(0.0, pX))
  } else {
    init <- control$init
  }

  # lower bound for optim
  lower <- c(rep(0.00001, 2*pW+1), rep(-Inf, pX))



  # call optimization
  if (is.null(control$cl)) {

  } else {
    stop("Parallelization not yet implemented.")
  }


  # tapering?
  taper <- if(is.null(control$taper)) {
    # without tapering
    NULL
  } else {
    # with tapering
    spam::cov.wend1(d, c(control$taper, 1, 0))
  }


  optim.output <- stats::optim(par = init,
                               fn = nLL,
                                 # arguments of nLL
                                 cov_func = cov.func,
                                 outer.W  = outer.W,
                                 y        = y,
                                 X        = X,
                                 W        = W,
                                 taper    = taper,
                               method = "L-BFGS-B",
                               lower = lower,
                               control = optim.control)

  # preparing output
  result <- list(optim.output = optim.output,
                 call.args = list(y = y,
                                  X = X,
                                  locs = locs,
                                  control = control,
                                  optim.control = optim.control, 
                                  W = W,
                                  ns = ns),
                 comp.args = list(outer.W = outer.W,
                                  lower   = lower,
                                  pW = pW,
                                  pX = pX,
                                  init = init))

  class(result) <- "SVC_mle"
  return(result)
}




#' Prediction of SVC (and response variable)
#'
#' @param object        output of \code{\link{SVC_mle}}
#' @param newlocs       matrix of dimension n' x 2. These are the new locations the SVCs are predicted for. If \code{NULL}, the locations from the \code{SVC_mle} (i.e. \code{locs}) are considered.
#' @param newX          optional matrix of dimension n' x pX. If provided, besides the predicted SVC, the function also returns the predicted response variable.
#' @param newW          optional matrix of dimension n' x pW.
#' @param backtransform logical. If standardization in function call of \code{SVC_mle} took place, backtransform results in prediciton. Does not cover transformation on model formula (e.g. log-transformation of response), only standardization of covariates.
#' @param ...           further arguments
#' @return returns a data frame of n' rows and with columns
#' \itemize{
#'   \item \code{SVC_1, ..., SVC_p}, i.e. the predicted SVC at locations \code{newlocs}
#'   \item \code{y.pred}, if \code{newX} is provided
#'   \item \code{loc_x, loc_y}, the locations of the predictions
#' }
#'
#' @seealso \code{\link{SVC_mle}}
#'
#' @import spam
#' @importFrom fields rdist
#' @importFrom stats dist sd
#' @export
predict.SVC_mle <- function(object, newlocs = NULL, newX = NULL, newW = NULL, backtransform = TRUE, ...) {

  hyper.par <- object$optim.output$par

  pW <- object$comp.args$pW
  pX <- object$comp.args$pX
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
    d <- as.matrix(stats::dist(object$call.args$locs))
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
                   p = pW,
                   cov_func = list(cov.func = cf, ns = object$call.args$ns),
                   outer.W = object$comp.args$outer.W)


  # cross-covariance beta' y
  cov_b_y <- Sigma_b_y(x = hyper.par,
                       cov.func = cf_dd,
                       W = as.matrix(object$call.args$W),
                       n.new = n.new)



  eff <- cov_b_y %*% spam::solve.spam(cov_y) %*%
    (object$call.args$y - object$call.args$X %*% hyper.par[2*pW + 1 + 1:pX])

  eff <- matrix(eff, ncol = pW)



  if (backtransform & ("scale.pars" %in% names(object$call.args$control))) {
    # effects back transformed
    eff.bt <- sapply(1:pW, function(j) {
      m.s <- object$call.args$control$scale.pars$scale.pars.W[, j]
      if (m.s[2] != 0) {
        return(eff[, j]/m.s[2])
      } else {
        return(eff[, j])
      }
    })
  }


  if (!is.null(newX) & !is.null(newW)) {

    stopifnot(pW == ncol(newW),
              pX == ncol(newX))

    if (backtransform & ("scale.pars" %in% names(object$call.args$control))) {
      newX <- sapply(1:pX, function(j) {
        m.s <- object$call.args$control$scale.pars$scale.pars.X[, j]
        if (m.s[2] != 0) {
          return((newX[, j]-m.s[1])/m.s[2])
        } else {
          return(newX[, j])
        }
      })

      newW <- sapply(1:pW, function(j) {
        m.s <- object$call.args$control$scale.pars$scale.pars.W[, j]
        if (m.s[2] != 0) {
          return((newW[, j]-m.s[1])/m.s[2])
        } else {
          return(newW[, j])
        }
      })

      correctionW <- sapply(1:pW, function(j) {
        m.s <- object$call.args$control$scale.pars$scale.pars.W[, j]
        if (m.s[2] != 0) {
          return(eff.bt[, j]*m.s[1])
        } else {
          return(0)
        }
      })

      correctionX <- sapply(1:pX, function(j) {
        m.s <- object$call.args$control$scale.pars$scale.pars.X[, j]
        if (m.s[2] != 0) {
          return((hyper.par[2*pW + 1 + 1:pX]*m.s[1])/m.s[2])
        } else {
          return(0)
        }
      })

      y.pred <- apply(as.matrix(newW) * eff.bt, 1, sum) +
        newX %*% hyper.par[2*pW + 1 + 1:pX] -
        sum(correctionX) -
        apply(correctionW, 1, sum)
    } else {
      y.pred <- apply(as.matrix(newW) * eff, 1, sum) + newX %*% hyper.par[2*pW + 1 + 1:pX]

    }


    out <- as.data.frame(cbind(eff, y.pred, newlocs))
    colnames(out) <- c(paste0("SVC_", 1:ncol(eff)), "y.pred", "loc_x", "loc_y")
  } else {
    out <- as.data.frame(cbind(eff, newlocs))
    colnames(out) <- c(paste0("SVC_", 1:ncol(eff)), "loc_x", "loc_y")
  }



  return(out)
}


