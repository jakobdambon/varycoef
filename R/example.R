
#' Sample Function for GP-based SVC Model for Given Locations
#'
#' @description Samples SVC data at given locations. The SVCs parameters and the
#' covariance function have to be provided. The sampled model matrix contains an
#' intercept as a first column and further covariates sampled from a standard
#' normal. The SVCs are sampled according to their given parametrization and at
#' respective observation locations. The error vector is sampled from a nugget
#' effect. Finally, the response vector is computed. Please note that the
#' function is not optimized for sampling large data sets.
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
#'    \code{locs} containing the original locations. The \code{true_pars}
#'    contains the data frame of covariance parameters that were used to
#'    sample the GP-based SVCs. The nugget variance has been added to the 
#'    original argument of the function with its respective variance, but 
#'    \code{NA} for \code{"mean"} and \code{"scale"}.
#'
#' @details The parameters of the model can be chosen such that we obtain data
#'    from a not full model, i.e., not all covariates are associated with a 
#'    fixed and a random effect. Using \code{var = 0} for instance yields a 
#'    constant beta coefficient for respective covariate. Note that in that 
#'    case the \code{scale} value is neglected.
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
#' SVCdata <- sample_SVCdata(
#'   df.pars = df.pars, nugget.sd = tau, locs = s, cov.name = "mat32"
#' )
#' @importFrom spam as.spam rmvnorm cov.exp cov.mat cov.sph cov.wend1 cov.wend2
#' @importFrom stats rnorm
#' @export
sample_SVCdata <- function(
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
  
  ## build SVC models depending on covariance function, i.e., Sigma_y
  D <- as.matrix(dist(locs, diag = TRUE, upper = TRUE))
  
  ## covariance functions
  cov_fun <- switch(
    match.arg(cov.name),
    "exp" = function(theta) {
      spam::cov.exp(D, theta = theta)
    },
    "sph" = function(theta) {
      spam::cov.sph(D, theta = theta)
    },
    "mat32" = function(theta) {
      spam::cov.mat(D, theta = c(theta, 3/2))
    },
    "mat52" = function(theta) {
      spam::cov.mat(D, theta = c(theta, 5/2))
    },
    "wend1" = function(theta) {
      spam::cov.wend1(D, theta = theta)
    },
    "wend2" = function(theta) {
      spam::cov.wend2(D, theta = theta)
    }
  )
  
  ## sample SVCs (including mean effect)
  beta <- apply(df.pars, 1, function(x) {
    if (x["var"] == 0) {
      rep(x["mean"], n)
    } else {
      spam::rmvnorm(
        n = 1,
        mu = rep(x["mean"], n), 
        Sigma = spam::as.spam(cov_fun(theta = x[c("scale", "var")]))
      )
    }
  })
  # nugget
  eps <- rnorm(n, sd = nugget.sd)
  
  # data
  X <- cbind(1, matrix(rnorm(n*(p-1)), ncol = p-1))
  y <- apply(beta*X, 1, sum) + eps
  
  list(
    y = y, X = X, beta = beta, eps = eps, locs = locs, 
    true_pars = rbind(
      df.pars, 
      data.frame(
        var = nugget.sd^2, 
        scale = NA, mean = NA
      )
    )
  )
}
