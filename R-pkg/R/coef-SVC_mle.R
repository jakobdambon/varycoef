
#' @title Extact Mean Effects
#'
#' @description Method to extract the mean effects from an \code{\link{SVC_mle}}
#' or \code{\link{SVC_selection}} object.
#'
#' @param object \code{\link{SVC_mle}} or \code{\link{SVC_selection}} object
#' @param ...    further arguments
#'
#' @return named vector with mean effects, i.e. \eqn{\mu} from
#' \code{\link[varycoef]{SVC_mle}}
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

#' @rdname coef.SVC_mle
#' @importFrom stats coef na.omit
#' @export
coef.SVC_selection <- function(object, ...) {
  mu <- as.numeric(tail(na.omit(object$PMLE_pars$mu.par), 1))

  X.vars <- colnames(object$obj.fun$args$X)

  names(mu) <- if (is.null(X.vars)) {
    paste0("Var", 1:length(mu))
  } else {
    X.vars
  }

  mu
}
