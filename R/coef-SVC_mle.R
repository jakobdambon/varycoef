
#' @title Extact Mean Effects
#'
#' @description Method to extract the mean effects from an \code{\link{SVC_mle}} object.
#'
#' @param object \code{\link{SVC_mle}} object
#' @param ...    further arguments
#'
#' @return named vector with mean effects, i.e. \eqn{\mu} from \code{\link[varycoef]{SVC_mle}}
#'
#' @author Jakob Dambon
#'
#' @importFrom stats coef
#' @export
coef.SVC_mle <- function(object, ...) {
  mu <- as.numeric(object$coefficients)

  X.vars <- colnames(object$data$X)

  names(mu) <- if (is.null(X.vars)) {
    paste0("Var", 1:length(mu))
  } else {
    X.vars
  }

  mu
}
