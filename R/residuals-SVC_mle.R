#' @title Extact Model Residuals
#'
#' @description Method to extract the residuals from an \code{\link{SVC_mle}} object. This is only possible if \code{save.fitted} was set to \code{TRUE} in the control of the function call
#'
#' @param object \code{\link{SVC_mle}} object
#' @param ...    further arguments
#'
#' @return numeric, residuals of model
#'
#' @author Jakob Dambon
#'
#' @importFrom stats residuals resid
#' @export
residuals.SVC_mle <- function(object, ...) {
  stopifnot(!is.null(object$residuals))
  return(object$residuals)
}

