#' Sample Function for SVCs
#'
#' @description Samples SVC on a regular grid. The SVC have all mean 0.
#'
#' @param m        integer. square root number of observations, in total the function will sample \eqn{m^2} locations on a regular grid.
#' @param p        integer. number of SVC
#' @param cov_pars data.frame including the covariance parameters of SVCs, using an exponential covariance function. The columns must have the names \code{"var"} and \code{"scale"}.
#' @param nugget   scalar. variance of the nugget / error term.
#' @param seed     integer. seed for sampling
#'
#' @return object of class \code{SpatialPointsDataFrame} (see \code{\link[sp]{SpatialPointsDataFrame-class}}) of the sampled SVC including the nugget.
#'
#' @examples 
#' # number of SVC
#' p <- 3
#' # sqrt of total number of observations
#' m <- 20
#' # covariance parameters
#' (pars <- data.frame(var = c(0.1, 0.2, 0.3), 
#'                     scale = c(0.3, 0.1, 0.2)))
#' nugget.var <- 0.05
#' 
#' # function to sample SVCs
#' sp.SVC <- fullSVC_reggrid(m = m, p = p, 
#'                           cov_pars = pars, 
#'                           nugget = nugget.var)
#' 
#' library(sp)
#' # visualization of sampled SVC
#' spplot(sp.SVC, colorkey = TRUE)
#'
#' @importFrom RandomFields RMnugget RFsimulate RMexp
#' @importFrom sp SpatialPointsDataFrame
#' @export
fullSVC_reggrid <- function(m, p, cov_pars, nugget, seed = 123) {

  # number of observations
  n <- m^2

  # regular grid locations
  locs <- expand.grid(x = seq(0, 1, length.out = m),
                      y = seq(0, 1, length.out = m))

  set.seed(seed)

  # SVC model
  model <- apply(cov_pars, 1, function(x) {
    RandomFields::RFsimulate(
      RandomFields::RMexp(x["var"], x["scale"]),
                          x = locs[, "x"], y = locs[, "y"])
  })

  model[[p+1]] <- RandomFields::RFsimulate(
    RandomFields::RMnugget(var = nugget),
                           x = locs[, "x"], y = locs[, "y"])
  sp.SVC <- Reduce(cbind, model)
  sp.SVC <- sp::SpatialPointsDataFrame(coords = sp.SVC@coords,
                                       data = sp.SVC@data)
  colnames(sp.SVC@data) <- c(paste0("SVC_", 1:p), "nugget")

  return(sp.SVC)
}
