
#' @title Extact Model Fitted Values
#'
#' @description Method to extract the fitted values from an \code{\link{SVC_mle}} object. This is only possible if \code{save.fitted} was set to \code{TRUE} in the control of the function call
#'
#' @param object \code{\link{SVC_mle}} object
#' @param ...    further arguments
#'
#' @return data frame, fitted values to given data, i.e. the SVC as well as the response and their locations
#'
#' @author Jakob Dambon
#'
#' @importFrom stats fitted
#' @export
fitted.SVC_mle <- function(object, ...) {
  stopifnot(!is.null(object$fitted))
  return(object$fitted)
}

