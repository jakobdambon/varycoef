#' Smoothly Clipped Absolute Deviation Penalty
#'
#' @description Penalty function proposed by Fan & Li (2001) \doi{10.1198/016214501753382273}.
#'
#' @param x      numeric.
#' @param lambda non-negative scalar, shrinkage parameter.
#' @param a      scalar larger than 2. Fan & Li (2001) suggest \eqn{a = 3.7}.
#'
#' @return penalty for values of \code{x}.
#' @export
#'
#' @author Jakob Dambon
#'
#' @references
#' Jianqing Fan & Runze Li (2001) Variable Selection via Nonconcave Penalized Likelihood and its Oracle Properties, Journal of the American Statistical Association, 96:456, 1348-1360, \doi{10.1198/016214501753382273}
#'
#' @examples
#'
#' SCAD(-5:5)
#' curve(SCAD, from = -5, to = 5)
SCAD <- function(x, lambda = 1, a = 3.7) {
  stopifnot(
    is.numeric(x) &
      lambda >= 0 & length(lambda) == 1 &
      a > 2 & length(a) == 1
  )

  ifelse(
    abs(x) <= lambda,
    lambda*abs(x),
    ifelse(
      abs(x) <= a*lambda,
      -(x^2 - 2*a*lambda*abs(x) + lambda^2)/(2*(a-1)),
      ((a+1)*lambda^2)/2))
}

#' Derivative of Smoothly Clipped Absolute Deviation Penalty
#'
#' @description derivative of \code{\link{SCAD}}, which is not differentiable at \code{x == 0}.
#'
#' @param x      numeric.
#' @param lambda non-negative scalar, shrinkage parameter.
#' @param a      scalar larger than 2. Fan & Li (2001) suggest \eqn{a = 3.7}.
#' @param d.side side of serivative at origin. Default value is \code{"both"}, returning NA for \code{x == 0}. If set to \code{"RHS"}, then returns RHS derivative, i.e., \eqn{\lambda}, and \eqn{-\lambda} with \code{"LHS"}.
#'
#' @return derivative of \code{SCAD(x)}, for \code{x == 0} return value is \code{NA}.
#' @export
#'
#' @author Jakob Dambon
#'
#' @examples
#'
#' d.SCAD(-5:5)
#' d.SCAD(-2:2, d.side = "LHS")
#' curve(d.SCAD, from = -5, to = 5)
d.SCAD <- function(x, lambda = 1, a = 3.7, d.side = "both") {
  stopifnot(
    is.numeric(x) &
      lambda >= 0 & length(lambda) == 1 &
      a > 2 & length(a) == 1
  )

  d.origin <- switch(
    d.side,
    "both" = NA,
    "RHS" = lambda,
    "LHS" = -lambda
  )


  ifelse(
     x == 0,
     d.origin, ifelse(
      abs(x) <= lambda,
      lambda*sign(x),
      ifelse(
        abs(x) <= a*lambda,
        sign(x)*(a*lambda-abs(x))/(a-1),
        0
      )
    )
  )
}



#' \eqn{L^q} Norm Penalty
#'
#' @description Penalty function using the \eqn{L^q} norm, i.e., \eqn{p_{\lambda}(x) = \lambda \| x\|^q}.
#'
#' @param x      numeric.
#' @param lambda non-negative scalar, shrinkage parameter.
#' @param q      non-negative scalar, norm parameter.
#'
#' @return penalty for values of \code{x}.
#' @export
#'
#' @author Jakob Dambon
#'
#' @examples
#'
#' Lq(-5:5)
#' curve(Lq(x, q = 2), from = -5, to = 5)
Lq <- function(x, lambda = 1, q = 1) {
  stopifnot(
    is.numeric(x) &
      lambda >= 0 & length(lambda) == 1 &
      q >= 0 & length(q) == 1
  )

  lambda*abs(x)^q
}

#' Derivative of \eqn{L^q} Norm Penalty
#'
#' @description derivative of \code{\link{Lq}}, which is not differentiable at \code{x == 0}.
#'
#' @param x      numeric.
#' @param lambda non-negative scalar, shrinkage parameter.
#' @param q      non-negative scalar, norm parameter..
#' @param d.side side of serivative at origin. Default value is \code{"both"}, returning NA for \code{x == 0}. If set to \code{"RHS"}, then returns RHS derivative, i.e., \eqn{\lambda}, and \eqn{-\lambda} with \code{"LHS"}.
#'
#' @return derivative of \code{Lq(x)}, i.e., \deqn{L^q(x) = q\lambda | x|^{q-1}}, for \code{x == 0} return value is \code{NA}.
#' @export
#'
#' @author Jakob Dambon
#'
#' @examples
#'
#' d.Lq(-5:5)
#' d.Lq(-2:2, d.side = "LHS")
#' curve(d.Lq, from = -5, to = 5)
d.Lq <- function(x, lambda = 1, q = 1, d.side = "both") {
  stopifnot(
    is.numeric(x) &
      lambda >= 0 & length(lambda) == 1 &
      q >= 0 & length(q) == 1
  )

  d.origin <- switch(
    d.side,
    "both" = NA,
    "RHS" = lambda,
    "LHS" = -lambda
  )

  ifelse(x == 0, d.origin, sign(x)*lambda*q*(abs(x)^(q-1)))
}
