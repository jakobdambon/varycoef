
#' @title Print Method for \code{SVC_mle}
#'
#' @description Method to print an \code{\link{SVC_mle}} object.
#'
#' @param x      \code{\link{SVC_mle}} object
#' @param digits  (\code{numeric})
#'    Number of digits to be plotted.
#' @param ...    further arguments
#'
#' @author Jakob Dambon
#'
#' @method print SVC_mle
#' @export
print.SVC_mle <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat(paste0("\nCall:\nSVC_mle with ",
             ncol(x$data$X),
             " fixed effect(s) and ",
             ncol(x$data$W),
             " SVC(s)\n"))
  cat(paste0("Using ", nobs(x), " observations at ",
             nlocs(x), " different locations.\n\n"))


  cat("Coefficients of fixed effects:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L,
                quote = FALSE)

  cat("\n\nCovaraiance parameters of the SVC(s):\n")
  print.default(format(cov_par(x), digits = digits), print.gap = 2L,
                quote = FALSE)

  cat("\n")
  invisible(x)

}
