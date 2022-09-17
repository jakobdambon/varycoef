
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

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")

  cat("Coefficients of fixed effects:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L,
                quote = FALSE)

  cat("\n\nCovaraiance parameters of the SVC(s):\n")
  print.default(format(cov_par(x), digits = digits), print.gap = 2L,
                quote = FALSE)

  cat("\n")
  invisible(x)

}
