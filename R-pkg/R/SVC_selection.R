BW_pen <- function(x, X, cov_func, outer.W, taper) {
  2*(eff_dof(x, X, cov_func, outer.W, taper) +
       2*length(outer.W) + 1)
}

VB_pen <- function(x, X, cov_func, outer.W, taper) {
  n <- nrow(X)
  p <- ncol(X)
  q <- length(outer.W)
  eff.dof <- eff_dof(x, X, cov_func, outer.W, taper)

  (2*n)/(n-p-2)*(eff.dof + 1 - (eff.dof - p)/(n-p))
}


#' @importFrom glmnet glmnet coef.glmnet
CD_mu <- function(
  theta.k,
  mle.par,
  obj.fun,
  lambda.mu,
  adaptive = FALSE
) {
  ## dimensions
  p <- ncol(obj.fun$args$X)
  q <- (length(theta.k)-1)/2

  ## transform from GLS to OLS
  # compute covariance matrix and Cholesky-decomp thereof
  C.mat <- Sigma_y(theta.k, obj.fun$args$cov_func, obj.fun$args$outer.W)
  R <- spam::chol(C.mat)
  R.t.inv <- solve(t(R))
  # transform
  y.tilde <- R.t.inv %*% obj.fun$args$y
  X.tilde <- R.t.inv %*% obj.fun$args$X

  # MLE / OLS estimate
  # mu.MLE <- coef(lm(y.tilde~X.tilde-1))

  ## run (adaptive) LASSO
  LASSO <- glmnet::glmnet(
    y = y.tilde,
    x = X.tilde,
    lambda = lambda.mu,
    alpha = 1,
    intercept = FALSE,
    penalty.factor = if (adaptive) {
      1/abs(mle.par)
    } else {
      rep(1, p)
    }
  )

  # without 0 for intercept
  as.numeric(coef(LASSO))[1+(1:p)]
}

CD_theta <- function(
  mu.k, theta.k, # theta.k only needed as initial value
  mle.par,
  obj.fun,
  lambda.theta,
  parallel.control,
  optim.args = list(),
  adaptive = FALSE
) {
  ## dimensions
  q <- (length(theta.k)-1)/2
  n <- length(obj.fun$args$y)

  # update profile log-lik. function using new mu.k
  obj.fun$args$mean.est <- mu.k

  fn <- function(x) {
    do.call(obj.fun$obj_fun, c(list(x = x), obj.fun$args))
  }

  # adaptive LASSO?
  if (adaptive) {
    # # MLE needed for adaptive penalties
    # MLE <- optimParallel::optimParallel(
    #   theta.k, fn = fn,
    #   lower = c(rep(c(1e-10, 0), q), 1e-10),
    #   parallel = parallel.control,
    #   control = list(
    #     parscale = ifelse(abs(theta.k) < 1e-9, 1, abs(theta.k))
    #   )
    # )
    # MLE.theta <- MLE$par


    # adaptive penalties
    lambda2 <- ifelse(
      abs(mle.par[2*(1:q)]) < 1e-9,
      1e99,
      lambda.theta/abs(mle.par[2*(1:q)])
    )
  } else {
    lambda2 <- lambda.theta
  }


  # objective function (f in paper)
  pl <- function(x) {
    fn(x) + 2*n*sum(lambda2*abs(x[2*(1:q)]))
  }

  # PMLE: use last known theta.k as initial value
  PMLE <- if (is.null(parallel.control)) {
    do.call(
      stats::optim,
      c(
        list(par = theta.k, fn = pl, method = "L-BFGS-B"),
        optim.args
      )
    )
  } else {
    do.call(
      optimParallel::optimParallel,
      c(
        list(
          par = theta.k,
          fn = pl,
          parallel = parallel.control
        ),
        optim.args
      )
    )
  }

  # return covariance parameters theta.k+1
  PMLE$par
}


PMLE_CD <- function(
  lambda,
  mle.par,
  obj.fun,
  parallel.control,
  optim.args = list(),
  adaptive = FALSE,
  return.par = FALSE,
  IC.type = c("BIC", "cAIC_BW", "cAIC_VB"),
  CD.conv = list(N = 20L, delta = 1e-6, logLik = TRUE)
) {
  ## dimensions
  n <- nrow(obj.fun$args$X)
  p <- ncol(obj.fun$args$X)
  q <- length(obj.fun$args$outer.W)

  ## initialize output matrix
  # covariance parameters
  c.par <- matrix(NA_real_, nrow = CD.conv$N + 1, ncol = 2*q+1)
  c.par[1, ] <- mle.par
  # mean parameter
  mu.par <- matrix(NA_real_, nrow = CD.conv$N + 1, ncol = p)
  I.C.mat <- solve(
    Sigma_y(mle.par, obj.fun$args$cov_func, obj.fun$args$outer.W)
  )
  B <- crossprod(obj.fun$args$X, I.C.mat)
  mu.par[1, ] <- solve(B %*% obj.fun$args$X) %*% B %*% obj.fun$args$y
  # log-likelihood
  loglik.CD <- rep(NA_real_, CD.conv$N + 1)

  # update mean parameter for log-likelihood function
  obj.fun$args$mean.est <- mu.par[1, ]
  # initialize log-likelihood function
  ll <- function(x) {
    (-1/2) * do.call(obj.fun$obj_fun, c(list(x = x), obj.fun$args))
  }
  loglik.CD[1] <- ll(c.par[1, ])

  ## cyclic coordinate descent
  for (k in 1:CD.conv$N) {

    # Step 1: Updating mu
    mu.par[k+1, ] <- CD_mu(
      theta.k = c.par[k, ],
      mle.par = mu.par[1, ],
      obj.fun = obj.fun,
      lambda.mu = lambda[1],
      adaptive = adaptive
    )

    # Step 2: Updating theta
    c.par[k+1, ] <- CD_theta(
      mu.k = mu.par[k+1, ],
      theta.k = c.par[k, ], # only used as an initial value for optimization
      mle.par = c.par[1, ],
      obj.fun = obj.fun,
      lambda.theta = lambda[2],
      parallel.control = parallel.control,
      optim.args = optim.args,
      adaptive = adaptive
    )

    ## compute new log-likelihood
    # update mean parameter for log-likelihood function
    obj.fun$args$mean.est <- mu.par[k + 1, ]
    # initialize log-likelihood function
    ll <- function(x) {
      (-1/2) * do.call(obj.fun$obj_fun, c(list(x = x), obj.fun$args))
    }
    loglik.CD[k + 1] <- ll(c.par[k+1, ])

    # check for convergence in theta parameters
    if (CD.conv$logLik) {
      # on the log likelihood
      if (abs(loglik.CD[k] - loglik.CD[k + 1])/
          abs(loglik.CD[k]) < CD.conv$delta) break
    } else {
      # on the parameters
      if (sum(abs(c(c.par[k, ], mu.par[k, ]) - c(c.par[k+1, ], mu.par[k+1, ])))/
          sum(abs(c(c.par[k, ], mu.par[k, ]))) < CD.conv$delta) break
    }
  }


  ## prepare output
  # update profile log-lik. function using last mu.k
  obj.fun$args$mean.est <- mu.par[k + 1, ]
  # note: neg2LL is the objective function of MLE is -2 times the log-lik.
  # We transform it back to the exact log-lik in the return call
  neg2LL <- function(x) {
    do.call(obj.fun$obj_fun, c(list(x = x), obj.fun$args))
  }

  # model-complexity (penalty)
  MC <- switch(
    match.arg(IC.type),
    cAIC_BW = BW_pen(
      c.par[k+1, ],
      obj.fun$args$X,
      obj.fun$args$cov_func,
      obj.fun$args$outer.W,
      obj.fun$args$taper
    ),
    cAIC_VB = VB_pen(
      c.par[k+1, ],
      obj.fun$args$X,
      obj.fun$args$cov_func,
      obj.fun$args$outer.W,
      obj.fun$args$taper
    ),
    BIC = {
      log(n)*sum(abs(c(mu.par[k+1, ], c.par[k+1, 2*(1:q)]) > 1e-10))
    }
  )
  # calculate final IC
  final.IC <- neg2LL(c.par[k+1, ]) + MC

  # return either all parameters of CD or
  # only IC value (needed for numeric optimization)
  if (return.par) {
    attr(final.IC, "IC.type") <- IC.type
    return(list(
      mu.par = mu.par,
      c.par = c.par,
      loglik.CD = loglik.CD,
      final.IC = final.IC
    ))
  } else {
    return(final.IC)
  }
}




#' @importFrom pbapply pbapply pboptions
IC_opt_grid <- function(IC.obj, r.lambda, n.lambda) {

  l.lambda <- seq(log(r.lambda[1]), log(r.lambda[2]), length.out = n.lambda)

  op <- pbapply::pboptions(type = "timer")
  IC_result <- pbapply::pbapply(
    expand.grid(lambda_mu = exp(l.lambda), lambda_sigma_sq = exp(l.lambda)),
    1,
    function(lambda)
      IC.obj(lambda = as.numeric(lambda))
  )
  pbapply::pboptions(op)

  out <- list(
    IC_grid = IC_result,
    l.lambda = l.lambda
  )
  class(out) <- c("SVC_pmle_grid", "SVC_pmle")
  return(out)
}



#' @importFrom pbapply pbapply pboptions
#' @importFrom ParamHelpers makeParamSet makeNumericVectorParam generateDesign
#' @importFrom smoof makeSingleObjectiveFunction
#' @importFrom mlr makeLearner
#' @importFrom mlrMBO makeMBOControl setMBOControlTermination mbo makeMBOInfillCritEI
#' @importFrom lhs maximinLHS
IC_opt_MBO <- function(
  IC.obj, r.lambda, n.init, n.iter,
  infill.crit = mlrMBO::makeMBOInfillCritEI()
) {

  par.set <- ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericVectorParam(
      "lambda",
      len = 2,
      lower = rep(r.lambda[1], 2),
      upper = rep(r.lambda[2], 2))
  )


  obj.fun <- smoof::makeSingleObjectiveFunction(
    fn = IC.obj,
    par.set = par.set,
    name = "IC"
  )

  design <- ParamHelpers::generateDesign(
    n = n.init,
    par.set = par.set,
    fun = lhs::maximinLHS
  )

  op <- pbapply::pboptions(type = "timer")
  design$y <- pbapply::pbapply(design, 1, obj.fun)
  pbapply::pboptions(op)

  surr.km <- mlr::makeLearner(
    "regr.km",
    predict.type = "se",
    covtype = "matern3_2"
  )

  control <- mlrMBO::makeMBOControl()
  control <- mlrMBO::setMBOControlTermination(control, iters = n.iter)
  control <- mlrMBO::setMBOControlInfill(
    control,
    crit = infill.crit
  )

  run <- mlrMBO::mbo(
    obj.fun,
    design = design,
    learner = surr.km,
    control = control,
    show.info = TRUE
  )

}


#' SVC Selection Parameters
#'
#' @description Function to set up control parameters for
#'    \code{\link{SVC_selection}}. The underlying Gaussian Process-based
#'    SVC model is defined in \code{\link{SVC_mle}}. \code{\link{SVC_selection}}
#'    then jointly selects fixed and random effects of the GP-based
#'    SVC model using a penalized maximum likelihood estimation (PMLE).
#'    In this function, one can set the parameters for the PMLE and
#'    its optimization procedures (Dambon et al., 2022).
#'
#' @param IC.type  (\code{character(1)}) \cr
#'    Select Information Criterion.
#' @param method   (\code{character(1)}) \cr
#'    Select optimization method for lambdas, i.e., shrinkage parameters.
#'    Either model-based optimization (MBO, Bischl et al., 2017 <arXiv:1703.03373>) or over grid.
#' @param r.lambda (\code{numeric(2)}) \cr
#'    Finite, positive, increasing lower and upper bounds for the shrinkage
#'    parameters.
#' @param n.lambda (\code{numeric(1)}) \cr
#'    If the grid method is selected, the positive integer number of values
#'    per side of the grid.
#' @param n.init   (\code{numeric(1)}) \cr
#'    If the MBO method is selected, the integer number of initial design
#'    points; at least three are required.
#' @param n.iter   (\code{numeric(1)}) \cr
#'    If the MBO method is selected, the positive integer number of
#'    surrogate-model iterations.
#' @param CD.conv  (\code{list(3)}) \cr
#'    Named list with positive integer \code{N}, finite positive \code{delta},
#'    and logical \code{logLik} elements. \code{N} is the maximum number of
#'    coordinate-descent iterations; \code{delta} is the relative convergence
#'    threshold; \code{logLik} selects convergence based on log-likelihood
#'    (\code{TRUE}) or parameter values (\code{FALSE}).
#' @param hessian  (\code{logical(1)}) \cr
#'    If \code{TRUE}, Hessian will be computed for final model.
#' @param adaptive (\code{logical(1)}) \cr
#'    If \code{TRUE}, adaptive LASSO is executed, i.e.,
#'    the shrinkage parameter is defined as \eqn{\lambda_j := \lambda / |\theta_j|}.
#' @param parallel (\code{list} or \code{NULL}) \cr
#'    Arguments for parallel optimization, as documented in
#'    \code{\link[optimParallel]{optimParallel}}. If \code{NULL}, serial
#'    optimization is used.
#' @param optim.args (\code{list}) \cr
#'    List of further arguments of \code{\link[optimParallel]{optimParallel}},
#'    such as the lower bounds.
#'
#' @export
#'
#' @examples
#' # Initializing parameters and switching logLik to FALSE
#' selection_control <- SVC_selection_control(
#'   CD.conv = list(N = 20L, delta = 1e-06, logLik = FALSE)
#' )
#' # or
#' selection_control <- SVC_selection_control()
#' selection_control$CD.conv$logLik <- FALSE
#'
#' @author Jakob Dambon
#'
#' @references Bischl, B., Richter, J., Bossek, J., Horn, D., Thomas, J.,
#'    Lang, M. (2017).
#'    \emph{mlrMBO: A Modular Framework for Model-Based Optimization of
#'    Expensive Black-Box Functions},
#'    ArXiv preprint \url{https://arxiv.org/abs/1703.03373}
#'
#'    Dambon, J. A., Sigrist, F., Furrer, R. (2022).
#'    \emph{Joint Variable Selection of both Fixed and Random Effects for
#'    Gaussian Process-based Spatially Varying Coefficient Models},
#'    International Journal of Geographical Information Science
#'    \doi{10.1080/13658816.2022.2097684}
#'
#'
#' @return A list of control parameters for SVC selection.
SVC_selection_control <- function(
  IC.type    = c("BIC", "cAIC_BW", "cAIC_VB"),
  method     = c("grid", "MBO"),
  r.lambda   = c(1e-10, 1e01),
  n.lambda   = 10L,
  n.init     = 10L,
  n.iter     = 10L,
  CD.conv    = list(N = 20L, delta = 1e-06, logLik = TRUE),
  hessian    = FALSE,
  adaptive   = FALSE,
  parallel   = NULL,
  optim.args = list()
) {

  IC.type <- match.arg(IC.type)
  method <- match.arg(method)

  if (!is.numeric(r.lambda) || length(r.lambda) != 2L ||
      any(!is.finite(r.lambda)) || any(r.lambda <= 0) ||
      r.lambda[1] >= r.lambda[2]) {
    stop("'r.lambda' must contain two finite, positive, increasing values.",
         call. = FALSE)
  }

  check_count <- function(x, name, minimum = 1L) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x != floor(x) || x < minimum) {
      stop(sprintf("'%s' must be a single integer greater than or equal to %d.",
                   name, minimum), call. = FALSE)
    }
  }
  check_count(n.lambda, "n.lambda")
  check_count(n.init, "n.init", minimum = if (method == "MBO") 3L else 1L)
  check_count(n.iter, "n.iter")

  if (!is.list(CD.conv) ||
      !all(c("N", "delta", "logLik") %in% names(CD.conv))) {
    stop("'CD.conv' must be a list with 'N', 'delta', and 'logLik' elements.",
         call. = FALSE)
  }
  check_count(CD.conv$N, "CD.conv$N")
  if (!is.numeric(CD.conv$delta) || length(CD.conv$delta) != 1L ||
      !is.finite(CD.conv$delta) || CD.conv$delta <= 0) {
    stop("'CD.conv$delta' must be a single finite positive number.",
         call. = FALSE)
  }
  if (!is.logical(CD.conv$logLik) || length(CD.conv$logLik) != 1L ||
      is.na(CD.conv$logLik)) {
    stop("'CD.conv$logLik' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(hessian) || length(hessian) != 1L || is.na(hessian)) {
    stop("'hessian' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(adaptive) || length(adaptive) != 1L || is.na(adaptive)) {
    stop("'adaptive' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(parallel) && !is.list(parallel)) {
    stop("'parallel' must be NULL or a list.", call. = FALSE)
  }
  if (!is.list(optim.args)) {
    stop("'optim.args' must be a list.", call. = FALSE)
  }


  list(
    IC.type    = IC.type,
    method     = method,
    r.lambda   = r.lambda,
    n.lambda   = n.lambda,
    n.init     = n.init,
    n.iter     = n.iter,
    CD.conv    = CD.conv,
    hessian    = hessian,
    adaptive   = adaptive,
    parallel   = parallel,
    optim.args = optim.args
  )
}


#' SVC Model Selection
#'
#' @description This function implements the variable selection for
#'   Gaussian process-based SVC models using a penalized maximum likelihood
#'   estimation (PMLE, Dambon et al., 2021, <arXiv:2101.01932>).
#'   It jointly selects the fixed and random effects of GP-based SVC models.
#'
#' @param obj.fun  (\code{SVC_obj_fun}) \cr
#'    Function of class \code{SVC_obj_fun}. This is the output of
#'    \code{\link{SVC_mle}} with the \code{\link{SVC_mle_control}} parameter
#'    \code{extract_fun} set to \code{TRUE}. This objective function comprises
#'    of the whole SVC model on which the selection should be applied.
#' @param mle.par   (\code{numeric(2*q+1)}) \cr
#'    Numeric vector with estimated covariance parameters of unpenalized MLE.
#' @param control   (\code{list} or \code{NULL}) \cr
#'    List of control parameters for variable selection. Output of
#'    \code{\link{SVC_selection_control}}. If \code{NULL} is given, the 
#'    default values of \code{\link{SVC_selection_control}} are used.
#' @param ...       Further arguments.
#'
#' @return Returns an object of class \code{SVC_selection}. It contains parameter estimates under PMLE and the optimization as well as choice of the shrinkage parameters.
#'
#' @author Jakob Dambon
#'
#' @references Dambon, J. A., Sigrist, F., Furrer, R. (2021).
#'    \emph{Joint Variable Selection of both Fixed and Random Effects for
#'    Gaussian Process-based Spatially Varying Coefficient Models},
#'    ArXiv Preprint \url{https://arxiv.org/abs/2101.01932}
#'
#' @export
#'
#'
#' @importFrom optimParallel optimParallel
SVC_selection <- function(
  obj.fun,
  mle.par,
  control = NULL,
  ...
) {

  if (!inherits(obj.fun, "SVC_obj_fun")) {
    stop("'obj.fun' must be an object of class 'SVC_obj_fun'; see ?SVC_mle.",
         call. = FALSE)
  }

  if (!is.list(obj.fun$args) ||
      !is.function(obj.fun$obj_fun) ||
      !is.matrix(obj.fun$args$X) ||
      !is.list(obj.fun$args$outer.W)) {
    stop("'obj.fun' is malformed: expected an objective function, a matrix 'X', ",
         "and a list 'outer.W'.", call. = FALSE)
  }

  n <- nrow(obj.fun$args$X)
  q <- length(obj.fun$args$outer.W)
  if (n < 1L || q < 1L) {
    stop("'obj.fun' must contain at least one observation and one random effect.",
         call. = FALSE)
  }

  if (!is.numeric(mle.par) || length(mle.par) != 2L * q + 1L ||
      any(!is.finite(mle.par))) {
    stop(sprintf(
      "'mle.par' must be a finite numeric vector of length %d.",
      2L * q + 1L
    ), call. = FALSE)
  }

  if (is.null(control)) {
    control <- SVC_selection_control()
  } else if (!is.list(control)) {
    stop("'control' must be NULL or a list returned by SVC_selection_control().",
         call. = FALSE)
  } else {
    control <- do.call(SVC_selection_control, control)
  }
  
  # IC black-box function
  IC.obj <- function(lambda)
    do.call(PMLE_CD, list(
      lambda = lambda,
      mle.par = mle.par,
      obj.fun = obj.fun,
      parallel.control = control$parallel,
      optim.args = control$optim.args,
      adaptive = control$adaptive,
      return.par = FALSE,
      IC.type = control$IC.type,
      CD.conv = control$CD.conv
    ))


  # start optimization
  PMLE_opt <- switch(
    control$method,
    "grid" = {IC_opt_grid(
      IC.obj = IC.obj,
      r.lambda = control$r.lambda,
      n.lambda = control$n.lambda
    )},
    "MBO"  = {
      stopifnot(control$n.init > 2)
      IC_opt_MBO(
        IC.obj = IC.obj,
        r.lambda = control$r.lambda,
        n.init = control$n.init,
        n.iter = control$n.iter
    )}
  )

  sel_lambda <- switch (
    control$method,
    "grid" = {
      l.grid <- expand.grid(
        PMLE_opt$l.lambda,
        PMLE_opt$l.lambda
      )

      exp(as.numeric(l.grid[which.min(PMLE_opt$IC_grid), ]))
    },
    "MBO"  = {
      as.numeric(PMLE_opt$x$lambda)
    }
  )


  PMLE <- function(lambda)
    do.call(PMLE_CD, list(
      lambda = lambda,
      mle.par = mle.par,
      obj.fun = obj.fun,
      parallel.control = control$parallel,
      optim.args = control$optim.args,
      adaptive = control$adaptive,
      return.par = TRUE,
      IC.type = control$IC.type,
      CD.conv = control$CD.conv
    ))

  PMLE_pars <- PMLE(sel_lambda)

  object <- list(
    PMLE_pars = PMLE_pars,
    PMLE_opt = PMLE_opt,
    lambda = sel_lambda,
    obj.fun = obj.fun,
    mle.par = mle.par
  )
  class(object) <- "SVC_selection"

  return(object)
}
