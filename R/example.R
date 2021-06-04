#' Sample Function for GP-based SVC Models on Regular Grid
#'
#' @description Samples SVC data on a regular quadratic (Cartesian) grid. The
#' SVCs have all mean 0 and an Matern covariance function is used.
#'
#' @param m   (\code{numeric(1)}) \cr
#'    Number of observations in one dimension, i.i, the square root number of
#'    total number of observation locations \eqn{n = m^2}.
#' @param p   (\code{numeric(1)}) \cr
#'    Number of SVCs.
#' @param cov_pars (\code{data.frame(p, 2)}) \cr
#'    Contains the covariance parameters of SVCs. The two columns must have the
#'    names \code{"nu"}, \code{"var"} and \code{"scale"}. These covariance
#'    parameters are then used for sampling the respective SVCs.
#' @param nugget  (\code{numeric(1)}) \cr
#'    Variance of the nugget / error term.
#' @param seed  (\code{numeric(1)}) \cr
#'    Seed set within the function for sampling.
#' @param given.locs (\code{NULL} or \code{data.frame(n, 2)}) \cr
#'    If \code{NULL}, the observations locations are sampled from a regular grid,
#'    Otherwise, the \code{data.frame} contains the observation locations.
#'    The data frame must have two columns of name \code{"x"} and \code{"y"}.
#'    The number of observations is then the number of rows \code{n}.
#'
#' @return (\code{data.frame}(m*m, p+3) \cr
#'    Data frame with \code{p+3} columns: first \code{p} are SVCs followed by a
#'    nugget and two columns for coordinates, i.e., locations. Number of
#'    observations is \code{m*m}.
#'
#' @examples
#' # number of SVC
#' p <- 3
#' # sqrt of total number of observations
#' m <- 20
#' # covariance parameters
#' (pars <- data.frame(
#'   nu = rep(0.5, p),
#'   var = c(0.1, 0.2, 0.3),
#'   scale = c(0.3, 0.1, 0.2)
#' ))
#' nugget.var <- 0.05
#'
#' # function to sample SVCs
#' sp.SVC <- fullSVC_reggrid(m = m, p = p,
#'                           cov_pars = pars,
#'                           nugget = nugget.var)
#' head(sp.SVC)
#'
#' @importFrom RandomFields RMnugget RFsimulate RMmatern
#' @export
fullSVC_reggrid <- function(m, p, cov_pars, nugget, seed = 123, given.locs = NULL) {

  if (is.null(given.locs)) {
    # number of observations
    n <- as.integer(m)^2

    # regular grid locations
    locs <- expand.grid(x = seq(0, 1, length.out = as.integer(m)),
                        y = seq(0, 1, length.out = as.integer(m)))
  } else {
    # take given locations
    locs <- given.locs
  }


  set.seed(seed)

  # SVC model
  model <- apply(cov_pars, 1, function(x) {
    RandomFields::RFsimulate(
      RandomFields::RMmatern(
        nu = x["nu"], var = x["var"], scale = x["scale"]),
      x = locs[, "x"], y = locs[, "y"])
  })

  model[[p+1]] <- RandomFields::RFsimulate(
    RandomFields::RMnugget(var = nugget),
                           x = locs[, "x"], y = locs[, "y"])
  SVCdata <- Reduce(cbind, model)
  SVCdata <- data.frame(data = SVCdata@data, coords = SVCdata@coords)
  colnames(SVCdata) <- c(paste0("SVC_", 1:p), "nugget", "locs1", "locs2")

  return(SVCdata)
}


#' Sample Function for GP-based SVC Model on Real Line
#'
#' @description Samples SVC data on a real line. The SVCs parameters and the
#' sample locations have to be provided. The SVCs are assumed to have an
#' Matern covariance function. The sampled model matrix contains an
#' intercept as a first column and further covariates sampled from a standard
#' normal. The SVCs are sampled according to their given parametrization and at
#' respective observation locations. The error vector sampled from a nugget
#' effect. Finally, the response vector is computed.
#'
#' @param df.pars (\code{data.frame(p, 3)}) \cr
#'    Contains the mean and covariance parameters of SVCs. The four columns
#'    must have the names \code{"mean"}, \code{"nu"}, \code{"var"}, and
#'    \code{"scale"}.
#' @param nugget.sd  (\code{numeric(1)}) \cr
#'    Standard deviation of the nugget / error term.
#' @param locs (\code{numeric(n)}) \cr
#'    The vector contains the observation locations and therefore defines the
#'    number of observations to be \code{n}.
#'
#' @return \code{list} \cr
#'    Returns a list with the response \code{y}, model matrix
#'    \code{X}, a matrix \code{beta} containing the sampled SVC at given
#'    locations, a vector \code{eps} containing the error, and a vector
#'    \code{locs} containing the original locations.
#'
#' @examples
#' set.seed(123)
#' # SVC parameters
#' (df.pars <- data.frame(
#'    nu = c(1.5, 1.5),
#'    var = c(2, 1),
#'    scale = c(3, 1),
#'    mean = c(1, 2)))
#' # nugget standard deviation
#' tau <- 0.5
#'
#' # sample locations
#' s <- sort(runif(500, min = 0, max = 10))
#' SVCdata <- fullSVC_line(
#'   df.pars = df.pars,
#'   nugget.sd = tau,
#'   locs = s
#' )
#'
#' @importFrom RandomFields RMtrend RFsimulate RMmatern
#' @importFrom stats rnorm
#' @export
fullSVC_line <- function(df.pars, nugget.sd, locs) {
  # dimensions
  n <- length(locs)
  p <- nrow(df.pars)
  # SVCs
  beta <- apply(df.pars, 1, function(x) {
    RandomFields::RFsimulate(
      RandomFields::RMmatern(nu = x["nu"], var = x["var"], scale = x["scale"]) +
        RandomFields::RMtrend(mean = x["mean"]),
      x = locs)@data$variable
  })
  # nugget
  eps <- rnorm(n, sd = nugget.sd)

  # data
  X <- cbind(1, matrix(rnorm(n*(p-1)), ncol = p-1))
  y <- apply(beta*X, 1, sum) + eps

  list(y = y, X = X, beta = beta, eps = eps, locs = locs)
}




#' Sample Function for GP-based SVC Model for Given Locations
#'
#' @description Samples SVC data at given locations. The SVCs parameters and the
#' covariance function have to be provided. The sampled model matrix contains an
#' intercept as a first column and further covariates sampled from a standard
#' normal. The SVCs are sampled according to their given parametrization and at
#' respective observation locations. The error vector sampled from a nugget
#' effect. Finally, the response vector is computed.
#'
#' @param df.pars (\code{data.frame(p, 3)}) \cr
#'    Contains the mean and covariance parameters of SVCs. The three columns
#'    must have the names \code{"mean"}, \code{"var"}, and \code{"scale"}.
#' @param nugget.sd  (\code{numeric(1)}) \cr
#'    Standard deviation of the nugget / error term.
#' @param cov.name (\code{character}(1)) \cr
#'    Character defining the covariance function, c.f. \code{\link{SVC_mle_control}}.
#' @param locs (\code{numeric(n)} or \code{matrix(n, d)}) \cr
#'    The numeric vector or matrix contains the observation locations and
#'    therefore defines the number of observations to be \code{n}. For a vector,
#'    we assume locations on the real line, i.e., \eqn{d=1}.
#'
#' @return \code{list} \cr
#'    Returns a list with the response \code{y}, model matrix
#'    \code{X}, a matrix \code{beta} containing the sampled SVC at given
#'    locations, a vector \code{eps} containing the error, and a matrix
#'    \code{locs} containing the original locations.
#'
#' @examples
#' set.seed(123)
#' # SVC parameters
#' (df.pars <- data.frame(
#'    var = c(2, 1),
#'    scale = c(3, 1),
#'    mean = c(1, 2)))
#' # nugget standard deviation
#' tau <- 0.5
#'
#' # sample locations
#' s <- sort(runif(500, min = 0, max = 10))
#' SVCdata <- sample_fullSVC(
#'   df.pars = df.pars, nugget.sd = tau, locs = s, cov.name = "mat32"
#' )
#'
#' @importFrom RandomFields RMtrend RFsimulate RMmatern RMgengneiting RMexp RMspheric
#' @importFrom stats rnorm
#' @export
sample_fullSVC <- function(
  df.pars, nugget.sd, locs,
  cov.name = c("exp", "sph", "mat32", "mat52", "wend1", "wend2")
) {
  # transform to matrix for further computations
  if (is.vector(locs)) {
    locs <- matrix(locs, ncol = 1)
  }
  # check covariance parameters and locations
  stopifnot(
    is.data.frame(df.pars),
    all(df.pars$var >= 0),
    all(df.pars$scale > 0),
    nugget.sd > 0,
    is.matrix(locs)
  )
  # dimensions
  d <- dim(locs)[2]
  n <- dim(locs)[1]
  p <- nrow(df.pars)
  # SVCs
  ## check Wendland and dimension
  if ((d>3)&(match.arg(cov.name) %in% c("wend1", "wend2"))) {
    warning("For Wendland covariance functions dimensions should be d<= 3.")
  }
  ## build SVC models depending on covariance function.
  RMmodels <- switch(
    match.arg(cov.name),
    "exp" = function(var, scale) {
      RandomFields::RMexp(var = var, scale = scale)
      },
    "sph" = function(var, scale) {
      RandomFields::RMspheric(var= var, scale = scale)
      },
    "mat32" = function(var, scale) {
      RandomFields::RMmatern(nu = 3/2, var = var, scale = scale)
      },
    "mat52" = function(var, scale) {
      RandomFields::RMmatern(nu = 5/2, var = var, scale = scale)
      },
    "wend1" = function(var, scale) {
      RandomFields::RMgengneiting (kappa = 1, var = var, scale = scale, mu = 3/2)
      },
    "wend2" = function(var, scale) {
      RandomFields::RMgengneiting (kappa = 2, var = var, scale = scale, mu = 3/2)
      }
  )
  ## sample SVCs (including mean effect)
  beta <- apply(df.pars, 1, function(x) {
    RandomFields::RFsimulate(
      RMmodels(var = x["var"], scale = x["scale"]) +
        RandomFields::RMtrend(mean = x["mean"]),
      x = locs)@data$variable
  })
  # nugget
  eps <- rnorm(n, sd = nugget.sd)

  # data
  X <- cbind(1, matrix(rnorm(n*(p-1)), ncol = p-1))
  y <- apply(beta*X, 1, sum) + eps

  list(y = y, X = X, beta = beta, eps = eps, locs = locs)
}


