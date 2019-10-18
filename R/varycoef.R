#' varycoef: Modeling Spatially Varying Coefficients
#'
#' This package offers functions to estimate and predict spatially varying coefficient (SVC) models. Briefly described, one generalizes a linear regression equation such that the coefficients are no longer constant, but have the possibility to vary spatially. This is enabled by modelling the coefficients by Gaussian random fields with either an exponential or spherical covariance function. The advantages of such SVC models are that they are usually quite easy to interpret, yet they offer a very highe level of flexibility.
#'
#' @section Estimation and Prediction:
#' The ensemble of the function \code{\link{SVC_mle}} and the method \code{predict} estimates the defined SVC model and gives predictions of the SVC as well as the response for some pre-defined locations. This concept should be rather familiar as it is the same for the classical regression (\code{\link{lm}}) or local polynomial regression (\code{\link{loess}}), to name a couple. As the name suggests, we are using a MLE approach in order to estimate the model and following the empirical best linear unbiased predictor to give location-specifc predictions. A detailed tutorial with examples is given in a vignette; call \code{vignette("example", package = "varycoef")}.
#'
#' @section Methods:
#' With the before mentioned \code{\link{SVC_mle}} function one gets an object of class \code{\link{SVC_mle}}. And like the method \code{predict} for predictions, there are several more methods in order to diagnose the model, see \code{methods(class = "SVC_mle")}.
#'
#' @examples
#' vignette("example", package = "varycoef")
#' methods(class = "SVC_mle")
#'
#' @docType package
#' @name varycoef
NULL
