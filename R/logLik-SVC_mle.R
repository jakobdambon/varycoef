

#' @title Extact the Likelihood
#'
#' @description Method to extract the computed (penalized) log (profile) Likelihood from an \code{\link{SVC_mle}} object.
#'
#' @param object \code{\link{SVC_mle}} object
#' @param ...    further arguments
#'
#' @return an object of class \code{logLik} with attributes
#' \itemize{
#'   \item \code{"penalized"}, logical, if the likelihood (\code{FALSE}) or some penalized likelihood (\code{TRUE}) was optimized.
#'   \item \code{"profileLik"}, logical, if the optimization was done using the profile likelihood  (\code{TRUE}) or not.
#'   \item \code{"nobs"}, integer of number of observations
#'   \item \code{"df"}, integer of how many parameters were estimated. \strong{Note}: This includes only the covariance parameters if the profile likelihood was used.
#' }
#'
#' @author Jakob Dambon
#'
#' @importFrom stats logLik
#' @export
logLik.SVC_mle <- function(object, ...) {

  profLik <- object$MLE$call.args$control$profileLik

  val <- (-1/2) * as.numeric(object$MLE$optim.output$value)
  attr(val, "penalized") <- (!is.null(object$MLE$call.args$control$pc.prior))
  attr(val, "profileLik") <- profLik
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- if (profLik) {
    length(cov_par(object))
  } else {
    length(cov_par(object)) + length(coef(object))
  }
  class(val) <- "logLik"
  val
}
