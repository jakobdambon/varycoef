#' @title Extact Covariance Parameters
#'
#' Function to extract the covariance parameters from an \code{\link{SVC_mle}} object.
#'
#' @param object \code{\link{SVC_mle}} object
#' @param ...    further arguments
#'
#' @return vector with covariance parameters and attributes what kind
#' \itemize{
#'   \item \code{"GRF"}, charachter, describing the GRF used, see \code{\link{SVC_mle_control}}.
#'   \item \code{"tapering"}, either \code{NULL} if no tapering is applied of the taper range.
#' }
#'
#' @author Jakob Dambon
#'
#'
#' @export
cov_par <- function(object, ...) {
  covpars <- as.numeric(object$cov.par)

  W.vars <- colnames(object$data$W)

  names(covpars) <- if (is.null(W.vars)) {
    c(paste0(rep(paste0("SVC", 1:((length(covpars)-1)/2)), each = 2),
             c(".range", ".var")), "nugget.var")
  } else {
    c(paste0(rep(W.vars, each = 2), c(".range", ".var")), "nugget.var")
  }

  attr(covpars, "GRF") <- object$MLE$call.args$control$cov.name
  attr(covpars, "tapering") <- object$MLE$call.args$control$tapering

  covpars
}









