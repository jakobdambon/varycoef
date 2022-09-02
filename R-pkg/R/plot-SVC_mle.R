

#' @title Plotting Residuals of \code{SVC_mle} model
#'
#' @description Method to plot the residuals from an \code{\link{SVC_mle}}
#' object. For this, \code{save.fitted} has to be \code{TRUE} in
#' \code{\link{SVC_mle_control}}.
#'
#' @param x          (\code{\link{SVC_mle}})
#' @param which      (\code{numeric}) \cr A numeric vector and subset of
#' \code{1:2} indicating which of the 2 plots should be plotted.
#' @param ...        further arguments
#'
#' @return a maximum 2 plots
#' \itemize{
#'   \item Tukey-Anscombe plot, i.e. residuals vs. fitted
#'   \item QQ-plot
#' }
#'
#' @author Jakob Dambon
#'
#' @seealso \code{\link[graphics]{legend}}  \link{SVC_mle}
#'
#' @examples
#' #' ## ---- toy example ----
#' ## sample data
#' # setting seed for reproducibility
#' set.seed(123)
#' m <- 7
#' # number of observations
#' n <- m*m
#' # number of SVC
#' p <- 3
#' # sample data
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p), ncol = p)
#' # locations on a regular m-by-m-grid
#' locs <- expand.grid(seq(0, 1, length.out = m),
#'                     seq(0, 1, length.out = m))
#'
#' ## preparing for maximum likelihood estimation (MLE)
#' # controls specific to MLE
#' control <- SVC_mle_control(
#'   # initial values of optimization
#'   init = rep(0.1, 2*p+1),
#'   # using profile likelihood
#'   profileLik = TRUE
#' )
#'
#' # controls specific to optimization procedure, see help(optim)
#' opt.control <- list(
#'   # number of iterations (set to one for demonstration sake)
#'   maxit = 1,
#'   # tracing information
#'   trace = 6
#' )
#'
#' ## starting MLE
#' fit <- SVC_mle(y = y, X = X, locs = locs,
#'                control = control,
#'                optim.control = opt.control)
#'
#' ## output: convergence code equal to 1, since maxit was only 1
#' summary(fit)
#'
#' ## plot residuals
#' # only QQ-plot
#' plot(fit, which = 2)
#'
#' # two plots next to each other
#' oldpar <- par(mfrow = c(1, 2))
#' plot(fit)
#' par(oldpar)
#'
#' @importFrom stats qqnorm qqline
#' @importFrom graphics plot abline legend
#' @method plot SVC_mle
#' @export
plot.SVC_mle <- function(x, which = 1:2, ...) {

  stopifnot(
    # residuals needed
    !is.null(x$residuals),
    # only two kinds of plots supported
    all(which %in% 1:2)
  )

  ## Tukey-Anscombe
  if (1 %in% which) {
    plot(fitted(x)$y.pred,
         residuals(x),
         main = "Tukey-Anscombe Plot",
         xlab = "fitted", ylab = "residuals",
         col = "grey")
    abline(h = 0)
  }

  ## QQ-plot
  if (2 %in% which) {
    qqnorm(residuals(x))
    qqline(residuals(x))
  }

  # ## spatial residuals
  # if (3 %in% which) {
  #   loc_x <- fitted(x)$loc_x
  #   loc_y <- fitted(x)$loc_y
  #   res <- residuals(x)
  #   cex.range <- range(sqrt(abs(res)))
  #
  #   plot(loc_x, loc_y,
  #        type = "p",
  #        main = "Spatial Residuals Plot",
  #        xlab = "x locations", ylab = "y locations",
  #        pch = 1, col = ifelse(res < 0, "blue", "orange"),
  #        cex = sqrt(abs(res)))
  #
  #   legend(legend.pos,
  #          legend = c("pos. residuals", "neg. residuals", "min", "max"),
  #          pch = c(19, 19, 1, 1),
  #          col = c("orange", "blue", "grey", "grey"),
  #          pt.cex = c(1, 1, cex.range))
  # }
}
