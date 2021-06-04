## -----------------------------------------------------------------------------
## In this script, one finds every function directly related to estimating
## and predicting SVC using our proposed MLE.
## -----------------------------------------------------------------------------



## ---- help function to give back correct covariance function ----
MLE.cov.func <- function(cov.name) {
  if (is.character(cov.name)) {
    cov.func <- switch(cov.name,
           "exp" = spam::cov.exp,
           "mat32" = function(h, theta) {
             spam::cov.mat(h, theta = c(theta, 3/2))},
           "mat52" = function(h, theta) {
             spam::cov.mat(h, theta = c(theta, 5/2))},
           "sph" = spam::cov.sph,
           "wend1" = spam::cov.wend1,
           "wend2" = spam::cov.wend2,
           stop("Cov.name argument not defined."))
  } else if (is.function(cov.name)) {
    cov.func <- cov.name
  } else {
    stop("Cov.name argument neither character, nor covariance function.")
  }
  return(cov.func)
}

## ---- help function to do MLE for SVC model ----
#' @importFrom stats coef lm median var
#' @importFrom optimParallel optimParallel
#' @importFrom parallel clusterExport clusterEvalQ
MLE_computation <- function(y, X, locs, W,
                            control,
                            optim.control) {
  ## -- set important dimensions ----
  # number random effects and fixed effects
  q <- dim(W)[2]
  p <- dim(X)[2]
  # indices of objective and covariance parameters
  id_obj <- if (control$profileLik) {
    (1:(2*q+1))
  } else {
    (1:(2*q+1+p))
  }
  id_cov <- (1:(2*q+1))

  # define distance matrix
  d <- do.call(
    own_dist,
    c(list(x = locs, taper = control$tapering), control$dist)
  )

  ## -- check and initialize optim vectors -----
  if (is.null(control$lower) | is.null(control$upper) | is.null(control$init)) {
    # median distances
    med_dist <- if (is.matrix(d)) {
      median(as.numeric(d))
    } else {
      options(spam.trivalues=TRUE)
      median(lower.tri(d@entries))
    }
    # variance of response
    y_var <- var(y)
    # fixed effects estimates by ordinary least squares (OLS)
    OLS_mu <- coef(lm(y~X-1))
  } else {
    med_dist <- y_var <- OLS_mu <- NULL
  }
  # liu - _L_ower _I_nit _U_pper
  liu <- init_bounds_optim(control, p, q, id_obj, med_dist, y_var, OLS_mu)

  ## -- define distance matrices, covariance functions, and taper matrix -----
  # get covariance function
  raw.cov.func <- MLE.cov.func(control$cov.name)

  # covariance function
  cov.func <- function(x) raw.cov.func(d, x)

  Rstruct <- NULL
  # tapering?
  if (is.null(control$tapering)) {
    taper <-NULL
    outer.W <- lapply(1:q, function(k) W[, k]%o%W[, k])
  } else {
    taper <- get_taper(control$cov.name, d, control$taper)
    outer.W <- lapply(1:q, function(k) {
      (W[, k]%o%W[, k]) * taper
    })

    options(spam.trivalues = TRUE, spam.cholsymmetrycheck = FALSE)

    Sigma1 <- Sigma_y(
      x = liu$init[id_cov],
      cov_func = cov.func,
      outer.W = outer.W,
      taper = taper
    )
    Rstruct <- spam::chol.spam(Sigma1)
  }

  ## -- pc priors -----
  # ordering: pcp = c(\rho_0, \alpha_\rho, \sigma_0, \alpha_\sigma)
  pcp.neg2dens <- if (is.null(control$pc.prior)) {
    NULL
  } else {
    pcp <- control$pc.prior

    lambda.r <- -log(pcp[2])*2*pcp[1]
    lambda.s <- -log(pcp[4])/pcp[3]

    # for Matérn GRF (-2 * log( pc prior dens))
    function(theta) {
      4*log(theta[1]) +
        lambda.r/theta[1]+2*lambda.s*sqrt(theta[2])
    }

  }


  # how to compute mu if optimization is over porfile likelihood
  # prepare for optimization by computing mean effect
  mu.estimate <- if (control$mean.est == "GLS") {
    NULL
  } else { # Ordinary Least Squares
    coef(lm(y~X-1))
  }

  # extract objective function
  if (control$extract_fun) {

    obj_fun <- function(x, ...)
      n2LL(x, ...)
    args <- list(
      cov_func = cov.func,
      outer.W  = outer.W,
      y        = y,
      X        = X,
      W        = W,
      mean.est = mu.estimate,
      taper    = taper,
      pc.dens  = pcp.neg2dens,
      Rstruct  = Rstruct,
      profile  = control$profileLik
    )

    return(list(
      obj_fun = obj_fun,
      args    = args
    ))
  }


  ## -- optimization -----
  if (is.null(control$parallel)) {
    # ... without parallelization
    optim.output <- stats::optim(
      par     = liu$init[id_obj],
      fn      = n2LL,
      # arguments of 2nLL
        cov_func = cov.func,
        outer.W  = outer.W,
        y        = y,
        X        = X,
        W        = W,
        mean.est = mu.estimate,
        taper    = taper,
        pc.dens  = pcp.neg2dens,
        Rstruct  = Rstruct,
        profile  = control$profileLik,
      method  = "L-BFGS-B",
      lower   = liu$lower[id_obj],
      upper   = liu$upper[id_obj],
      hessian = control$hessian,
      control = optim.control)
  } else {
    # ... with parallelization
    parallel::clusterEvalQ(
      cl = control$parallel$cl,
      {
        library(spam)
        library(varycoef)
      }
    )

    parallel::clusterExport(
      cl = control$parallel$cl,
      varlist = ls(),
      envir = environment()
    )

    optim.output <- optimParallel::optimParallel(
      par     = liu$init[id_obj],
      fn      = n2LL,
      # arguments of 2nLL
        cov_func = cov.func,
        outer.W  = outer.W,
        y        = y,
        X        = X,
        W        = W,
        mean.est = mu.estimate,
        taper    = taper,
        pc.dens  = pcp.neg2dens,
        Rstruct  = Rstruct,
        profile  = control$profileLik,
      lower   = liu$lower[id_obj],
      upper   = liu$upper[id_obj],
      hessian = control$hessian,
      control = optim.control,
      parallel = control$parallel
    )
  }

  ## -- Estimates and Standard errors -----
  # compute covariance matrices
  Sigma_final <- Sigma_y(
    optim.output$par[id_cov], cov.func, outer.W, taper = taper
  )

  par_SE <- prep_par_output(
    optim.output$par, Sigma_final, Rstruct, control$profileLik, X, y,
    optim.output$hessian, q
  )

  # effective degrees of freedom
  edof <- eff_dof(
    cov.par  = par_SE$RE$est,
    cov_func = cov.func,
    outer.W  = outer.W,
    X        = X,
    taper    = taper
  )

  # preparing output
  return(
    list(
      optim.output = optim.output,
      call.args = list(
        y = as.numeric(y),
        X = as.matrix(X),
        locs = as.matrix(locs),
        control = control,
        optim.control = optim.control,
        W = W
      ),
      comp.args = list(
        liu         = liu,
        edof        = edof,
        Sigma_final = Sigma_final,
        par_SE      = par_SE
      )
    )
  )
}

## ---- help function to compute fitted values after MLE ----
fitted_computation <- function(SVC_obj, y, X, W, locs) {
  class(SVC_obj) <- "SVC_mle"


  predict.SVC_mle(SVC_obj, newlocs = locs, newX = X, newW = W)

}

## ---- help function to construct SVC_mle object ----
create_SVC_mle <- function(ML_estimate, y, X, W, locs, control) {

  q <- dim(W)[2]

  # extract covariance parameters and coefficients for methods
  cov.par <- ML_estimate$comp.args$par_SE$RE$est
  mu <- ML_estimate$comp.args$par_SE$FE$est

  # non zero parameters, i.e., means or variances
  df <- sum(abs(c(mu, cov.par[2*(1:q)])) > 1e-10)

  SVC_obj <- list(
    MLE          = ML_estimate,
    coefficients = mu,
    cov.par      = cov.par,
    # (effective) degrees of freedom
    df = list(
      df = as.integer(df),
      edof = ML_estimate$comp.args$edof),
    fitted = NULL,
    residuals = NULL,
    data = list(y = y, X = X, W = W, locs = locs)
  )


  if (control$save.fitted) {
    # compute fitted values (i.e. EBLUP = empirical BLUP)
    pred <- fitted_computation(SVC_obj, y, X, W, locs)

    SVC_obj$fitted = pred
    SVC_obj$residuals = y-pred$y.pred
  }

  return(SVC_obj)
}


#' @title Set Parameters for \code{SVC_mle}
#'
#' @description Function to set up control parameters for \code{\link{SVC_mle}}.
#' In the following, we assume the GP-based SVC model to have \eqn{q} GPs which
#' model the SVCs and \eqn{p} fixed effects.
#'
#' @param cov.name  (\code{character(1)}) \cr
#'    Name of the covariance function of the GPs. Currently, the following are
#'    implemented: \code{"exp"} for the exponential, \code{"sph"} for
#'    spherical, \code{"mat32"} and \code{"mat52"} for Matern class covariance
#'    functions with smoothness 3/2 or 5/2, as well as \code{"wend1"} and
#'    \code{"wend2"} for Wendland class covariance functions with kappa 1 or 2.
#' @param tapering  (\code{NULL} or \code{numeric(1)}) \cr
#'    If \code{NULL}, no tapering is applied. If a scalar is given, covariance
#'    tapering with this taper range is applied, for all Gaussian processes
#'    modeling the SVC. Only defined for Matern class covariance functions,
#'    i.e., set \code{cov.name} either to \code{"exp"}, \code{"mat32"}, or
#'    \code{"mat52"}.
#' @param parallel  (\code{NULL} or \code{list}) \cr
#'    If \code{NULL}, no parallelization is applied. If cluster has been
#'    established, define arguments for parallelization with a list, see
#'    documentation of \code{\link[optimParallel]{optimParallel}}. See Examples.
#' @param init  (\code{NULL} or \code{numeric(2q+1+p*as.numeric(profileLik))}) \cr
#'    Initial values for optimization procedure. If \code{NULL} is given, an
#'    initial vector is calculated (see Details). Otherwise, the vector is
#'    assumed to consist of q-times (alternating) range and variance,
#'    the nugget variance and if \code{profileLik = TRUE} p mean effects.
#' @param lower  (\code{NULL} or \code{numeric(2q+1+p*as.numeric(profileLik))}) \cr
#'    Lower bound for \code{init} in \code{optim}. Default \code{NULL} calculates
#'    the lower bounds (see Details).
#' @param upper  (\code{NULL} or \code{numeric(2q+1+p*as.numeric(profileLik))}) \cr
#'    Upper bound for \code{init} in \code{optim}. Default \code{NULL} calculates
#'    the upper bounds (see Details).
#' @param save.fitted (\code{logical(1)}) \cr
#'    If \code{TRUE}, calculates the fitted values and residuals after MLE and
#'    stores them. This is necessary to call \code{\link{residuals}} and
#'    \code{\link{fitted}} methods afterwards.
#' @param profileLik  (\code{logical(1)}) \cr
#'    If \code{TRUE}, MLE is done over profile Likelihood of covariance
#'    parameters.
#' @param mean.est  (\code{character(1)}) \cr
#'    If \code{profileLik = TRUE}, the means have to be estimated seperately for
#'    each step. \code{"GLS"} uses the generalized least square estimate while
#'    \code{"OLS"} uses the ordinary least squares estimate.
#' @param pc.prior  (\code{NULL} or \code{numeric(4)}) \cr
#'    If numeric vector is given, penalized complexity priors are applied. The
#'    order is \eqn{\rho_0, \alpha_\rho, \sigma_0, \alpha_\sigma} to give some
#'    prior believes for the range and the standard deviation of GPs, such that
#'    \eqn{P(\rho < \rho_0) = \alpha_\rho, P(\sigma > \sigma_0) = \alpha_\sigma}.
#'    This regulates the optimization process. Currently, only supported for
#'    GPs with of Matérn class covariance functions. Based on the idea by
#'    Fulgstad et al. (2018) \doi{10.1080/01621459.2017.1415907}.
#' @param extract_fun (\code{logical(1)}) \cr
#'    If \code{TRUE}, the function call of \code{\link{SVC_mle}} stops before
#'    the MLE and gives back the objective function of the MLE as well as all
#'    used arguments. If \code{FALSE}, regular MLE is conducted.
#' @param hessian  (\code{logical(1)}) \cr
#'    If \code{TRUE}, Hessian matrix is computed, see \link[stats]{optim}. This
#'    required to give the standard errors for covariance parameters and to do
#'    a Wald test on the variances, see \code{\link{summary.SVC_mle}}.
#' @param dist     (\code{list}) \cr
#'    List containing the arguments of \link[stats]{dist} or
#'    \link[spam]{nearest.dist}. This controls
#'    the method of how the distances and therefore dependency structures are
#'    calculated. The default gives Euclidean distances in a \eqn{d}-dimensional
#'    space. Further editable arguments are \code{p, miles, R}, see respective
#'    help files of \link[stats]{dist} or \link[spam]{nearest.dist}.
#' @param ...     Further Arguments yet to be implemented
#'
#' @details If not provided, the initial values as well as the lower and upper
#'    bounds are calculated given the provided data. In particular, we require
#'    the median distance between observations, the variance of the response and,
#'    the ordinary least square (OLS) estimates, see \code{\link{init_bounds_optim}}.
#'
#'    The argument \code{extract_fun} is useful, when one wants to modify
#'    the objective function. Further, when trying to parallelize the
#'    optimization, it is useful to check whether a single evaluation of the
#'    objective function takes longer than 0.05 seconds to evaluate,
#'    cf. Gerber and Furrer (2019) \doi{10.32614/RJ-2019-030}. Platform specific
#'    issues can be sorted out by the user by setting up their own optimization.
#'
#' @return A list with which \code{\link{SVC_mle}} can be controlled.
#' @seealso \code{\link{SVC_mle}}
#'
#' @examples
#' control <- SVC_mle_control(init = rep(0.3, 10))
#' # or
#' control <- SVC_mle_control()
#' control$init <- rep(0.3, 10)
#'
#' \donttest{
#' # Code for setting up parallel computing
#' require(parallel)
#' # exchange number of nodes (1) for detectCores()-1 or appropriate number
#' cl <- makeCluster(1, setup_strategy = "sequential")
#' clusterEvalQ(
#'   cl = cl,
#'   {
#'     library(spam)
#'     library(varycoef)
#'   })
#' # use this list for parallel argument in SVC_mle_control
#' parallel.control <- list(cl = cl, forward = TRUE, loginfo = TRUE)
#' # SVC_mle goes here ...
#' # DO NOT FORGET TO STOP THE CLUSTER!
#' stopCluster(cl); rm(cl)
#' }
#' @author Jakob Dambon
#'
#' @export
SVC_mle_control <- function(...) UseMethod("SVC_mle_control")


#' @rdname SVC_mle_control
#' @export
SVC_mle_control.default <- function(
  cov.name = c("exp", "sph", "mat32", "mat52", "wend1", "wend2"),
  tapering = NULL,
  parallel = NULL,
  init = NULL,
  lower = NULL,
  upper = NULL,
  save.fitted = TRUE,
  profileLik = FALSE,
  mean.est = c("GLS", "OLS"),
  pc.prior = NULL,
  extract_fun = FALSE,
  hessian = TRUE,
  dist = list(method = "euclidean"),
  ...
) {
  stopifnot(
    is.null(tapering) | (tapering>=0),
    is.logical(save.fitted),
    is.logical(profileLik),
    is.logical(extract_fun),
    is.logical(hessian)
  )

  # if (!is.null(tapering) &
  #     !(match.arg(cov.name) %in% c("sph", "wend1", "wend2"))) {
  #   stop("Covariance tapering only defined for Matern class covariance functions.")
  # }

  list(
    cov.name = match.arg(cov.name),
    tapering = tapering,
    parallel = parallel,
    init = init,
    lower = lower,
    upper = upper,
    save.fitted = save.fitted,
    profileLik = profileLik,
    mean.est = match.arg(mean.est),
    pc.prior = pc.prior,
    extract_fun = extract_fun,
    hessian = hessian,
    dist = dist,
    ...
  )
}

#' @param object  (\code{SVC_mle}) \cr
#'   The function then extracts the control settings from the function call
#'   used to compute in the given \code{SVC_mle} object.
#'
#' @rdname SVC_mle_control
#' @export
SVC_mle_control.SVC_mle <- function(object, ...) {
  object$MLE$call.args$control
}





###############################
## SVC MLE functions ##########
###############################


#' @title MLE of SVC model
#'
#' @description Conducts a maximum likelihood estimation (MLE) for a Gaussian
#'   process-based SVC model as described in Dambon et al. (2021)
#'   \doi{10.1016/j.spasta.2020.100470}. More specifially, the model is
#'   defined as:
#'
#' \deqn{y(s) = X \mu + W \eta (s) + \epsilon(s)}
#'
#' where:
#' \itemize{
#'   \item \eqn{y} is the response (vector of length \eqn{n})
#'   \item \eqn{X} is the data matrix for the fixed effects covariates. The
#'   dimensions are \eqn{n} times \eqn{p}. This leads to \eqn{p} fixed effects.
#'   \item \eqn{\mu} is the vector containing the fixed effects
#'   \item W is the data matrix for the SVCs modeled by GPs. The dimensions are
#'   \eqn{n} times \eqn{q}. This lead to \eqn{q} SVCs in the model.
#'   \item \eqn{\eta} are the SVCs represented by a GP.
#'   \item \eqn{\epsilon} is the nugget effect
#' }
#'
#' The MLE is an numeric optimization that runs \code{\link[stats]{optim}} or
#' (if parallelized) \code{\link[optimParallel]{optimParallel}}.
#'
#' @param y  (\code{numeric(n)}) \cr
#'    Response vector.
#' @param X  (\code{matrix(n, p)}) \cr
#'    Design matrix. Intercept has to be added manually.
#' @param locs  (\code{matrix(n, d)}) \cr
#'    Locations in a \eqn{d}-dimensional space. May contain multiple
#'    observations at single location.
#' @param W  (\code{NULL} or \code{matrix(n, q)}) \cr
#'    If \code{NULL}, the same matrix as provided in \code{X} is used. This
#'    fits a full SVC model, i.e., each covariate effect is modeled with a mean
#'    and an SVC. In this case we have \eqn{p = q}. If optional matrix \code{W}
#'    is provided, SVCs are only modeled for covariates within matrix \code{W}.
#' @param control  (\code{list}) \cr
#'    Control paramaters given by \code{\link{SVC_mle_control}}.
#' @param optim.control  (\code{list}) \cr
#'    Control arguments for optimization function, see Details in
#'    \code{\link{optim}}.
#' @param ...            further arguments
#'
#' @return Object of class \code{SVC_mle} if \code{control$extract_fun = FALSE},
#' meaning that a MLE has been conducted. Otherwise, if \code{control$extract_fun = TRUE},
#' the function returns a list with two entries:
#' \itemize{
#'    \item \code{obj_fun}: the objective function used in the optimization
#'    \item \code{args}: the arguments to evaluate the objective function.
#' }
#' For further detials, see description of \code{\link{SVC_mle_control}}.
#'
#' @references Dambon, J. A., Sigrist, F., Furrer, R. (2021)
#'    \emph{Maximum likelihood estimation of spatially varying coefficient
#'    models for large data with an application to real estate price prediction},
#'    Spatial Statistics \doi{10.1016/j.spasta.2020.100470}
#' @author Jakob Dambon
#'
#' @seealso \code{\link{predict.SVC_mle}}
#'
#' @examples
#' ## ---- toy example ----
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
#'   # lower bound
#'   lower = rep(1e-6, 2*p+1),
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
#' class(fit)
#'
#' ## output: convergence code equal to 1, since maxit was only 1
#' summary(fit)
#'
#' ## extract the optimization arguments, including objective function
#' control$extract_fun <- TRUE
#' opt <- SVC_mle(y = y, X = X, locs = locs,
#'                control = control)
#'
#' # objective function and its arguments of optimization
#' class(opt$obj_fun)
#' class(opt$args)
#'
#' # single evaluation with initial value
#' do.call(opt$obj_fun,
#'         c(list(x = control$init), opt$args))
#'
#' \donttest{
#' ## ---- real data example ----
#' require(sp)
#' ## get data set
#' data("meuse", package = "sp")
#'
#' # construct data matrix and response, scale locations
#' y <- log(meuse$cadmium)
#' X <- model.matrix(~1+dist+lime+elev, data = meuse)
#' locs <- as.matrix(meuse[, 1:2])/1000
#'
#'
#' ## starting MLE
#' # the next call takes a couple of seconds
#' fit <- SVC_mle(y = y, X = X, locs = locs,
#'                # has 4 fixed effects, but only 3 random effects (SVC)
#'                # elev is missing in SVC
#'                W = X[, 1:3],
#'                control = SVC_mle_control(
#'                  # inital values for 3 SVC
#'                  # 7 = (3 * 2 covariance parameters + nugget)
#'                  init = c(rep(c(0.4, 0.2), 3), 0.2),
#'                  profileLik = TRUE
#'                ))
#'
#' ## summary and residual output
#' summary(fit)
#' plot(fit)
#'
#' ## predict
#' # new locations
#' newlocs <- expand.grid(
#'   x = seq(min(locs[, 1]), max(locs[, 1]), length.out = 30),
#'   y = seq(min(locs[, 2]), max(locs[, 2]), length.out = 30))
#' # predict SVC for new locations
#' SVC <- predict(fit, newlocs = as.matrix(newlocs))
#' # visualization
#' sp.SVC <- SVC
#' coordinates(sp.SVC) <- ~loc_1+loc_2
#' spplot(sp.SVC, colorkey = TRUE)
#' }
#' @import spam
#' @importFrom stats dist optim
#' @importFrom optimParallel optimParallel
#' @export
SVC_mle <- function(...) UseMethod("SVC_mle")



#' @rdname SVC_mle
#' @export
SVC_mle.default <- function(y, X, locs, W = NULL,
                            control = NULL,
                            optim.control = list(), ...) {

  # check if W is given arguments
  if (is.null(W)) {W <- X}

  # call SVC_mle with default control settings if non are provided
  if (is.null(control)) {
    control <- SVC_mle_control()
  }

  # Start ML Estimation using optim
  ML_estimate <- MLE_computation(y = y,
                                 X = X,
                                 locs = locs,
                                 W = W,
                                 control = control,
                                 optim.control = optim.control)

  if (is.function(ML_estimate$obj_fun)) {
    # extract objective function
    object <- ML_estimate
    class(object) <- "SVC_obj_fun"
    return(object)
  } else {
    # after optimization
    object <- create_SVC_mle(ML_estimate, y, X, W, locs, control)
    object$call <- match.call()
    class(object) <- "SVC_mle"
    return(object)
  }

}

# formula call

#' @param formula Formula describing the fixed effects in SVC model. The response, i.e. LHS of the formula, is not allowed to have functions such as \code{sqrt()} or \code{log()}.
#' @param data data frame containing the observations
#' @param RE_formula Formula describing the random effects in SVC model. Only RHS is considered. If \code{NULL}, the same RHS of argument \code{formula} for fixed effects is used.
#' @importFrom stats model.matrix
#'
#' @rdname SVC_mle
#' @export
SVC_mle.formula <- function(formula, data, RE_formula = NULL,
                            locs, control, optim.control = list(), ...) {


  X <- as.matrix(model.matrix(formula, data = data))
  W <- if (is.null(RE_formula)) {X} else {
    as.matrix(model.matrix(RE_formula, data = data))
  }
  y <- as.numeric(data[, all.vars(formula)[1]])

  SVC_mle.default(y = y,
                  X = X,
                  locs = locs,
                  W = W,
                  control = control,
                  optim.control = optim.control)
}

