#' @title Conditional Akaike's and Bayesian Information Criteria
#'
#' @name IC.SVC_mle
#' @aliases AIC.SVC_mle
#'
#' @description Methods to calculate information criteria for
#' \code{\link{SVC_mle}} objects. Currently, two are supported: the conditional
#' Akaike's Information Criteria \eqn{cAIC = -2*log-likelihood + 2*(edof + df)}
#' and the Bayesian Information Criteria \eqn{BIC = -2*log-likelihood + log(n) * npar}.
#' Note that the Akaike's Information Criteria is of the corrected form, that
#' is: \eqn{edof} is the effective degrees of freedom which is derived as the
#' trace of the hat matrices and df is the degree of freedoms with respect to
#' mean parameters.
#'
#' @param object \code{\link{SVC_mle}} object
#' @param ...    further arguments
#'
#' @return numeric, value of information criteria
#'
#' @author Jakob Dambon
#' @importFrom stats BIC
#' @export
BIC.SVC_mle <- function(object, ...) {
  -2*logLik(object) +
    log(nobs(object)) *
    object$df$df
}


#' @rdname IC.SVC_mle
#' @param conditional string. If \code{conditional = "BW"}, the
#' conditional AIC is calculated.
#'
#' @importFrom stats AIC
#' @export
AIC.SVC_mle <- function(object, conditional = "BW", ...) {
  -2*logLik(object) +
    2*(object$df$edof + object$df$df)
}

