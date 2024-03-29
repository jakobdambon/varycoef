

#' Prediction of SVCs (and response variable)
#'
#' @param object  (\code{SVC_mle}) \cr
#'    Model obtained from \code{\link{SVC_mle}} function call.
#' @param newlocs  (\code{NULL} or \code{matrix(n.new, 2)}) \cr
#'    If \code{NULL}, then function uses observed locations of model to estimate
#'    SVCs. Otherwise, these are the new locations the SVCs are predicted for.
#' @param newX  (\code{NULL} or \code{matrix(n.new, q)}) \cr
#'    If provided (together with \code{newW}), the function also returns the
#'    predicted response variable.
#' @param newW  (\code{NULL} or \code{matrix(n.new, p)}) \cr
#'    If provided (together with \code{newX}), the function also returns the
#'    predicted response variable.
#' @param newdata (\code{NULL} or \code{data.frame(n.new, p)}) \cr
#'    This argument can be used, when the \code{SVC_mle} function has been called
#'    with an formula, see examples.
#' @param compute.y.var  (\code{logical(1)}) \cr
#'    If \code{TRUE} and the response is being estimated, the predictive
#'    variance of each estimate will be computed.
#' @param ...           further arguments
#'
#' @return The function returns a data frame of \code{n.new} rows and with
#' columns
#' \itemize{
#'   \item \code{SVC_1, ..., SVC_p}: the predicted SVC at locations \code{newlocs}.
#'   \item \code{y.pred}, if \code{newX} and \code{newW} are provided
#'   \item \code{y.var}, if \code{newX} and \code{newW} are provided and
#'   \code{compute.y.var} is set to \code{TRUE}.
#'   \item \code{loc_x, loc_y}, the locations of the predictions
#' }
#'
#' @seealso \code{\link{SVC_mle}}
#'
#' @author Jakob Dambon
#' @references Dambon, J. A., Sigrist, F., Furrer, R. (2021)
#'    \emph{Maximum likelihood estimation of spatially varying coefficient
#'    models for large data with an application to real estate price prediction},
#'    Spatial Statistics \doi{10.1016/j.spasta.2020.100470}
#'
#' @examples
#' ## ---- toy example ----
#' ## We use the sampled, i.e., one dimensional SVCs
#' str(SVCdata)
#' # sub-sample data to have feasible run time for example
#' set.seed(123)
#' id <- sample(length(SVCdata$locs), 50)
#' 
#' ## SVC_mle call with matrix arguments
#' fit_mat <- with(SVCdata, SVC_mle(
#'   y[id], X[id, ], locs[id], 
#'   control = SVC_mle_control(profileLik = TRUE, cov.name = "mat32")))
#' 
#' ## SVC_mle call with formula
#' df <- with(SVCdata, data.frame(y = y[id], X = X[id, -1]))
#' fit_form <- SVC_mle(
#'   y ~ X, data = df, locs = SVCdata$locs[id], 
#'   control = SVC_mle_control(profileLik = TRUE, cov.name = "mat32")
#' )
#'
#' ## prediction
#' 
#' # predicting SVCs
#' predict(fit_mat, newlocs = 1:2)
#' predict(fit_form, newlocs = 1:2)
#' 
#' # predicting SVCs and response providing new covariates
#' predict(
#'   fit_mat, 
#'   newX = matrix(c(1, 1, 3, 4), ncol = 2), 
#'   newW = matrix(c(1, 1, 3, 4), ncol = 2), 
#'   newlocs = 1:2
#' )
#' predict(fit_form, newdata = data.frame(X = 3:4), newlocs = 1:2)
#' 
#' @import spam
#' @importFrom stats sd model.matrix
#' @export
predict.SVC_mle <- function(
  object,
  newlocs = NULL,
  newX = NULL,
  newW = NULL,
  newdata = NULL,
  compute.y.var = FALSE,
  ...
) {
  # extract parameters
  mu <- coef(object)
  cov.par <- cov_par(object)


  q <- dim(as.matrix(object$MLE$call.args$W))[2]
  p <- dim(as.matrix(object$MLE$call.args$X))[2]
  n <- length(object$MLE$call.args$y)
  locs <- object$MLE$call.args$locs
  tapering <- object$MLE$call.args$control$tapering
  dist_args <- object$MLE$call.args$control$dist

  # define distance matrices
  d <- do.call(
    own_dist,
    c(list(x = locs, taper = tapering), dist_args)
  )

  # if no new locations are given, predict for training data
  if (is.null(newlocs)) {
    newlocs <- object$MLE$call.args$locs
    d_cross <- d
    n.new <- n
  } else {
    newlocs <- as.matrix(newlocs)
    n.new <- nrow(newlocs)
    d_cross <- do.call(
      own_dist,
      c(list(x = newlocs, y = locs, taper = tapering), dist_args)
    )
  }

  # covariance function (not tapered)
  raw.cf <- MLE.cov.func(object$MLE$call.args$control$cov.name)

  if (is.null(object$MLE$call.args$control$taper)) {
    taper <- NULL

    # cross-covariance (newlocs and locs)
    cf_cross <- function(x) raw.cf(d_cross, x)

  } else {
    taper <- get_taper(
      object$MLE$call.args$control$cov.name, d, tapering
    )

    taper_cross <- get_taper(
      object$MLE$call.args$control$cov.name, d_cross, tapering
    )

    # cross-covariance (newlocs and locs)
    cf_cross <- function(x) raw.cf(d_cross, x)*taper_cross
  }

  # covariance y
  cf <- function(x) raw.cf(d, x)
  cov_y <- object$MLE$comp.args$Sigma_final


  # cross-covariance beta' y
  cov_b_y <- Sigma_b_y(
    x = cov.par,
    cov.func = cf_cross,
    W = as.matrix(object$MLE$call.args$W),
    n.new = n.new
  )

  eff <- cov_b_y %*%
    solve(cov_y, object$MLE$call.args$y - object$MLE$call.args$X %*% mu)
  eff <- matrix(eff, ncol = q)

  # if newdata is given and formula is present in SVC_mle object, extract
  # newX and newW (and overwrite provided ones)
  if (!is.null(newdata)) {
    if (!is.null(object$formula)) {
      if (!is.null(newX)) {
        warning("Formula and 'newdata' provided: 'newX' argument was overwritten!")
      }
      if (!is.null(newW)) {
        warning("Formula and 'newdata' provided: 'newW' argument was overwritten!")
      }
      # create covariates
      # drop response from fromula
      formula <- drop_response(object$formula)
      RE_formula <- drop_response(object$RE_formula)
      
      newX <- as.matrix(stats::model.matrix(formula, data = newdata))
      newW <- as.matrix(stats::model.matrix(RE_formula, data = newdata))
    } else {
      warning("Data provided bu object has not been trained by a formula.\n
              Cannot compute fixed and random effect covariates.")
    }
  }

  if (!is.null(newX) & !is.null(newW)) {
    # Do dimensions for training and prediction data match?
    stopifnot(q == ncol(newW), p == ncol(newX))
    y.pred <- apply(newW * eff, 1, sum) + newX %*% mu

    # computation of standard deviation fro each observation.
    if (compute.y.var) {
      # Have to compute
      #
      # var.y = Sigma_ynew - Sigma_ynew_y Sigma_y^-1 Sigma_y_ynew
      #
      # Sigma_ynew   = A
      # Sigma_ynew_y = B
      # Sigma_y      = C
      # Sigma_y_ynew = D = t(C)


      # Part B:
      cov_ynew_y <- Sigma_y_y(
        cov.par,
        cov.func = cf_cross,
        X = object$MLE$call.args$W,
        newX = newW
      )

      # Part A:
      d_new <- if (n.new == 1) {
        as.matrix(0)
      } else {
        do.call(
          own_dist,
          c(list(x = newlocs, taper = tapering), dist_args)
        )
      }

      if (is.null(tapering)) {
        outer.newW <- lapply(1:q, function(k) {
          (newW[, k]%o%newW[, k]) })
        taper_new <- NULL
      } else {
        taper_new <- get_taper(
          object$MLE$call.args$control$cov.name, d_new, tapering
        )

        outer.newW <- lapply(1:q, function(k) {
          (newW[, k]%o%newW[, k]) * taper_new
        })
      }

      # cross-covariance (newlocs and locs)
      cf_new <- function(x) raw.cf(d_new, x)

      cov_ynew <- Sigma_y(
        cov.par,
        cf_new,
        outer.W = outer.newW,
        taper = taper_new
      )

      # Part C: already calculated with cov_y

      # Computation of variance of y
      var.y <- diag(cov_ynew) - diag(cov_ynew_y %*% solve(cov_y, t(cov_ynew_y)))

      # form out put
      out <- as.data.frame(cbind(eff, y.pred, var.y, newlocs))
      colnames(out) <- c(paste0("SVC_", 1:ncol(eff)), "y.pred", "y.var", paste0("loc_", 1:ncol(newlocs)))
    } else {
      out <- as.data.frame(cbind(eff, y.pred, newlocs))
      colnames(out) <- c(paste0("SVC_", 1:ncol(eff)), "y.pred", paste0("loc_", 1:ncol(newlocs)))
    }


  } else {

    if (compute.y.var)
      warning("Please provide 'newX' and 'newW' to predict y and its variance.")

    out <- as.data.frame(cbind(eff, newlocs))
    colnames(out) <- c(paste0("SVC_", 1:ncol(eff)), paste0("loc_", 1:ncol(newlocs)))
  }
  # two ensure that predict calls with formula and matrix are identical
  row.names(out) <- as.character(row.names(out))
  return(out)
}

