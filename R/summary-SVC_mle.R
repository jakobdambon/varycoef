
#' @title Summary Method for \code{SVC_mle}
#'
#' @description Method to construct a \code{summary.SVC_mle} object out of a
#' \code{\link{SVC_mle}} object.
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
#' @importFrom stats pchisq pnorm
#' @method summary SVC_mle
#' @export
summary.SVC_mle <- function(object, ...) {

  stopifnot(!is.null(object$residuals))

  p <- dim(as.matrix(object$data$X))[2]
  q <- dim(as.matrix(object$data$W))[2]

  se_RE <- object$MLE$comp.args$par_SE$RE$SE
  se_FE <- object$MLE$comp.args$par_SE$FE$SE

  covpars <- cbind(
    Estimate = cov_par(object),
    `Std. Error` = se_RE,
    `W value` = (cov_par(object)/se_RE)^2,
    `Pr(>W)` = pchisq((cov_par(object)/se_RE)^2, df = 1, lower.tail = FALSE)
  )
  # do not test range and nugget variance
  covpars[-(2*(1:q)), 3:4] <- NA

  ans <- list(
    pX = p,
    pW = q,
    nobs = nobs(object),
    nlocs = nlocs(object),
    resids = resid(object),
    y.mean.resid = object$data$y-mean(object$data$y),
    coefs = cbind(
      Estimate = coef(object),
      `Std. Error` = se_FE,
      `Z value` = (coef(object)/se_FE),
      `Pr(>|Z|)` = 2*pnorm(abs(coef(object)/se_FE), lower.tail = FALSE)
    ),
    covpars = covpars,
    cov_fun = switch(
      attr(cov_par(object), "cov_fun"),
      "exp" = "exponential",
      "mat32" = "Matern (nu = 3/2)",
      "mat52" = "Matern (nu = 5/2)",
      "sph" = "spherical",
      "wend1" = "Wendland (kappa = 1)",
      "wend2" = "Wendland (kappa = 2)"),
    optim.out = object$MLE$optim.output,
    logLik = logLik(object),
    taper = object$MLE$call.args$control$tapering,
    BIC = as.numeric(BIC(object))
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
#' @importFrom stats printCoefmat sd
#' @method print summary.SVC_mle
#' @export
print.summary.SVC_mle <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  cat(paste0("\nCall:\nSVC_mle with ",
             x$pX,
             " fixed effect(s) and ",
             x$pW,
             " SVC(s)\n"))
  cat(paste0("using ", x$nobs, " observations at ",
             x$nlocs, " different locations / coordinates.\n\n"))

  cat("Residuals:\n")
  print.default(format(summary(x$resids)[-4], digits = digits), print.gap = 2L,
                quote = FALSE)
  cat(paste0("\nResidual standard error: ",
             formatC(sd(x$resids),
                     digits = digits),
             "\nMultiple R-squared: ",
             formatC(x$r.squared,
                     digits = digits),
             ", BIC: ", formatC(x$BIC,
                                digits = digits), "\n"))

  cat("\n\nCoefficients of fixed effect(s):\n")
  stats::printCoefmat(x$coefs, digits = digits,
                      signif.stars = getOption("show.signif.stars"),
                      na.print = "NA", ...)
  # print.default(format(x$coefs, digits = digits), print.gap = 2L,
  #               quote = FALSE)

  # covpar <- as.data.frame(matrix(x$covpars[-(2*x$pW+1)],
  #                                ncol = 2, byrow = TRUE))
  # colnames(covpar) <- c("range", "variance")
  # rownames(covpar) <- substr(names(x$covpars)[2*(1:x$pW)],
  #                            1, nchar(names(x$covpars)[2*(1:x$pW)])-4)

  cat("\n\nCovariance parameters of the SVC(s):\n")
  stats::printCoefmat(x$covpar, digits = digits,
                      signif.stars = getOption("show.signif.stars"),
                      na.print = "NA", ...)

  cat(paste0("\nThe covariance parameters were estimated using \n",
             x$cov_fun, " covariance functions.\n"))
  if(is.null(x$taper)) {
    cat("No covariance tapering applied.\n")
  } else {
    cat(paste0("Covariance tapering range set to: ",
               formatC(x$taper, digits = digits), "\n"))
  }


  cat(paste0(
    "\n\nMLE:\nThe MLE terminated after ",
    x$optim.out$counts["function"],
    " function evaluations with convergence code ",
    x$optim.out$convergence,
    "\n(0 meaning that the optimization was succesful).\n"
  ))
  cat(paste0(
    "The final", if (attr(x$logLik, "penalized")) {" regularized"} else {""} ,
    if (attr(x$logLik, "profileLik")) {" profile"} else {""},
    " log likelihood value is " ,
    formatC(x$logLik,digits = digits), ".\n"
  ))
  cat("\n")
  invisible(x)

}
