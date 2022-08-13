#' @title Extact Covariance Parameters
#'
#' @description Function to extract the covariance parameters from an 
#' \code{\link{SVC_mle}} or \code{\link{SVC_selection}}object.
#'
#' @param object \code{\link{SVC_mle}} or \code{\link{SVC_selection}} object
#' @param ...    further arguments
#'
#' @return vector with covariance parameters with the following attributes:
#' \itemize{
#'   \item \code{"GRF"}, charachter, describing the covariance function used for
#'   the GP, see \code{\link{SVC_mle_control}}.
#'   \item \code{"tapering"}, either \code{NULL} if no tapering is applied of
#'   the taper range.
#' }
#'
#' @author Jakob Dambon
#'
#'
#' @export
cov_par <- function(...) UseMethod("cov_par")

#' @rdname cov_par
#' @export
cov_par.SVC_mle <- function(object, ...) {
  covpars <- as.numeric(object$cov.par)

  W.vars <- colnames(object$data$W)

  names(covpars) <- if (is.null(W.vars)) {
    c(paste0(rep(paste0("SVC", 1:((length(covpars)-1)/2)), each = 2),
             c(".range", ".var")), "nugget.var")
  } else {
    c(paste0(rep(W.vars, each = 2), c(".range", ".var")), "nugget.var")
  }

  attr(covpars, "cov_fun") <- object$MLE$call.args$control$cov.name
  attr(covpars, "tapering") <- object$MLE$call.args$control$tapering

  covpars
}


#' @rdname cov_par
#' @importFrom stats na.omit
#' @export
cov_par.SVC_selection <- function(object, ...) {
  covpars <- as.numeric(tail(na.omit(object$PMLE_pars$c.par), 1))

  W.vars <- colnames(object$obj.fun$args$W)

  names(covpars) <- if (is.null(W.vars)) {
    c(paste0(rep(paste0("SVC", 1:((length(covpars)-1)/2)), each = 2),
             c(".range", ".var")), "nugget.var")
  } else {
    c(paste0(rep(W.vars, each = 2), c(".range", ".var")), "nugget.var")
  }

  attr(covpars, "cov_fun") <- attr(object$mle.par, "cov_fun")
  attr(covpars, "tapering") <- object$obj.fun$args$taper

  covpars
}






