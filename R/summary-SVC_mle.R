
#' @title Summary Method for \code{SVC_mle}
#'
#' @description Method to construct a \code{summary.SVC_mle} object out of a \code{\link{SVC_mle}} object.
#'
#' @param object \code{\link{SVC_mle}} object
#' @param ...    further arguments
#'
#' @return object of class \code{summary.SVC_mle} with summarized values of the MLE.
#'
#' @author Jakob Dambon
#'
#' @seealso \code{\link{SVC_mle}}
#'
#' @method summary SVC_mle
#' @export
summary.SVC_mle <- function(object, ...) {

  stopifnot(!is.null(object$residuals))

  ans <- list(
    pX = ncol(object$data$X),
    pW = ncol(object$data$W),
    nobs = nobs(object),
    nlocs = nlocs(object),
    resids = resid(object),
    y.mean.resid = object$data$y-mean(object$data$y),
    coefs = coef(object),
    covpars = cov_par(object),
    optim.out = object$MLE$optim.output,
    logLik = logLik(object),
    taper = object$MLE$call.args$control$tapering
  )

  ans$r.squared <- 1 - sum(ans$resids^2)/sum(ans$y.mean.resid^2)
  class(ans) <- "summary.SVC_mle"
  ans

}


#' @title Printing Method for \code{summary.SVC_mle}
#'
#' @param x      \code{\link{summary.SVC_mle}}
#' @param digits the number of significant digits to use when printing.
#' @param ...    further arguments
#'
#'
#' @return The printed output of the summary in the console.
#' @seealso \link{summary.SVC_mle} \link{SVC_mle}
#'
#' @method print summary.SVC_mle
#' @export
print.summary.SVC_mle <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(paste0("\nCall:\nSVC_mle with ",
             x$pX,
             " fixed effect(s) and ",
             x$pW,
             " SVC(s)\n"))
  cat(paste0("using ", x$nobs, " observations at ",
             x$nlocs, " different locations.\n\n"))

  cat("Residuals:\n")
  print.default(format(summary(x$resids)[-4], digits = digits), print.gap = 2L,
                quote = FALSE)
  cat(paste0("\nResidual standard error (nugget effect): ",
             formatC(sqrt(x$covpars[2*x$pW+1]),
                     digits = digits),
             "\nMultiple R-squared: ",
             formatC(x$r.squared,
                     digits = digits), "\n"))

  cat("\n\nCoefficients of fixed effects:\n")
  print.default(format(x$coefs, digits = digits), print.gap = 2L,
                quote = FALSE)

  covpar <- as.data.frame(matrix(x$covpars[-(2*x$pW+1)],
                                 ncol = 2, byrow = TRUE))
  colnames(covpar) <- c("range", "variance")
  rownames(covpar) <- substr(names(x$covpars)[2*(1:x$pW)],
                             1, nchar(names(x$covpars)[2*(1:x$pW)])-4)

  cat("\n\nCovariance parameters of the SVC(s):\n")
  print(covpar, digits = digits)

  cat(paste0("\nThe covariance parameters were estimated for the GRFs using\n",
             attr(x$covpars, "GRF"), ". covariance functions. "))
  if(is.null(x$taper)) {
    cat("No covariance tapering applied.\n")
  } else {
    cat(paste0("Covariance tapering range set to: ",
               formatC(x$taper, digits = digits), "\n"))
  }


  cat(paste0("\n\nMLE:\nThe MLE terminated after ",
             x$optim.out$counts["function"], " function evaluations with convergence code ",
             x$optim.out$convergence, "\n(0 meaning that the optimization was succesful).\n"))
  cat(paste0("The final", if (attr(x$logLik, "penalized")) {" regularized"} else {""} ,
             if (attr(x$logLik, "profileLik")) {" profile"} else {""},
             " likelihood value for ", attr(x$logLik, "df"), " parameters is " ,
             formatC(x$logLik,digits = digits), ".\n"))
  cat("\n")
  invisible(x)

}
