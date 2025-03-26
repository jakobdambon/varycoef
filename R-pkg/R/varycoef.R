#' varycoef: Modeling Spatially Varying Coefficients
#'
#' This package offers functions to estimate and predict Gaussian process-based
#' spatially varying coefficient (SVC) models. Briefly described, one
#' generalizes a linear regression equation such that the coefficients are no
#' longer constant, but have the possibility to vary spatially. This is enabled
#' by modeling the coefficients using Gaussian processes with (currently) either
#' an exponential or spherical covariance function. The advantages of such SVC
#' models are that they are usually quite easy to interpret, yet they offer a
#' very high level of flexibility.
#'
#'
#' @section Estimation and Prediction:
#' The ensemble of the function \code{\link{SVC_mle}} and the method
#' \code{predict} estimates the defined SVC model and gives predictions of the
#' SVC as well as the response for some pre-defined locations. This concept
#' should be rather familiar as it is the same for the classical regression
#' (\code{\link{lm}}) or local polynomial regression (\code{\link{loess}}),
#' to name a couple. As the name suggests, we are using a \emph{maximum
#' likelihood estimation} (MLE) approach in order to estimate the model. The
#' predictor is obtained by the empirical best linear unbiased predictor.
#' to give location-specific predictions. A detailed tutorial with examples is
#' given in a vignette; call \code{vignette("example", package = "varycoef")}.
#' We also refer to the original article Dambon et al. (2021) which lays the
#' methodological foundation of this package.
#'
#'
#' With the before mentioned \code{\link{SVC_mle}} function one gets an object
#' of class \code{\link{SVC_mle}}. And like the method \code{predict} for
#' predictions, there are several more methods in order to diagnose the model,
#' see \code{methods(class = "SVC_mle")}.
#'
#' @section Variable Selection:
#' As of version 0.3.0 of \code{varycoef}, a joint variable selection of both
#' fixed and random effect of the Gaussian process-based SVC model is
#' implemented. It uses a \emph{penalized maximum likelihood estimation} (PMLE)
#' which is implemented via a gradient descent. The estimation of the shrinkage
#' parameter is available using a \emph{model-based optimization} (MBO). Here,
#' we use the framework by Bischl et al. (2017). The methodological foundation
#' of the PMLE is described in Dambon et al. (2022).
#'
#' @examples
#' vignette("manual", package = "varycoef")
#' methods(class = "SVC_mle")
#'
#' @author Jakob Dambon
#'
#' @references Bischl, B., Richter, J., Bossek, J., Horn, D., Thomas, J.,
#'    Lang, M. (2017). \emph{mlrMBO: A Modular Framework for Model-Based
#'    Optimization of Expensive Black-Box Functions},
#'    ArXiv preprint \url{https://arxiv.org/abs/1703.03373}
#'
#'    Dambon, J. A., Sigrist, F., Furrer, R. (2021).
#'    \emph{Maximum likelihood estimation of spatially varying coefficient
#'    models for large data with an application to real estate price prediction},
#'    Spatial Statistics 41 100470 \doi{10.1016/j.spasta.2020.100470}
#'
#'    Dambon, J. A., Sigrist, F., Furrer, R. (2022).
#'    \emph{Joint Variable Selection of both Fixed and Random Effects for
#'    Gaussian Process-based Spatially Varying Coefficient Models},
#'    International Journal of Geographical Information Science
#'    \doi{10.1080/13658816.2022.2097684}
#'
#' @docType package
#' @name varycoef
#' @aliases varycoef-package
#' @rdname varycoef
NULL
