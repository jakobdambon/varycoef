#' @title Extract Number of Observations
#'
#' @description Method to extract the number of observations used in MLE for an \code{\link{SVC_mle}} object.
#'
#' @param object \code{\link{SVC_mle}} object
#' @param ...    further arguments
#'
#' @return an integer of number of observations
#'
#' @author Jakob Dambon
#'
#' @importFrom stats nobs
#' @export
nobs.SVC_mle <- function(object, ...) {
  length(object$data$y)
}

