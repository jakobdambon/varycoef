## ----setup, include=FALSE------------------------------------------------
rm(list = ls())
knitr::opts_chunk$set(echo = TRUE)
library(varycoef)

## ----define and sample SVC, warning=FALSE, fig.width=7, fig.height=7, message=FALSE----
# number of SVC
p <- 3

(pars <- data.frame(mu = rep(0, p), 
                    var = c(0.1, 0.2, 0.3), 
                    scale = c(0.3, 0.1, 0.2)))
nugget.var <- 0.05

# sqrt of total number of observations

m <- 20


library(RandomFields)
library(sp)

sp.SVC <- varycoef:::fullSVC_reggrid(m = m, p = p, 
                                     pars = list(pars = pars, 
                                                 nugget.var = nugget.var))

spplot(sp.SVC, colorkey = TRUE)

## ----sample covariates---------------------------------------------------
n <- m^2

X <- matrix(c(rep(1, n), rnorm((p-1)*n)), ncol = p)
head(X)

## ----compute response----------------------------------------------------
y <- apply(X * as.matrix(sp.SVC@data[, 1:p]), 1, sum) + sp.SVC@data[, p+1]

## ----contorl parameters--------------------------------------------------
control <- SVC_mle.control()
str(control)

## ----mean initials-------------------------------------------------------
(mu.init <- coef(lm(y~.-1, data = data.frame(y = y, X = X))))

## ----variance initials---------------------------------------------------
(var.init <- sigma(lm(y~.-1, data = data.frame(y = y, X = X)))^2)

## ----scale initials------------------------------------------------------
scale.init <- 0.4

## ----joint initials------------------------------------------------------
init <- c(rep(c(scale.init, var.init), p), # GRFs scales and variances
          var.init,                        # nugget variance
          mu.init)                         # means

## ----overwrite initials--------------------------------------------------
# default
control$init

# overwrite
control$init <- init

# create new
control <- SVC_mle.control(init = init)

## ----SVC MLE-------------------------------------------------------------
fit <- SVC_mle(y = y, X = X, locs = coordinates(sp.SVC), control = control)

class(fit)

# comparison of estimated and true parameters
rbind(fit$optim.output$par, 
      c(pars[, "scale"], pars[, "var"], nugget.var, pars[, "mu"]))

## ----make predictions----------------------------------------------------
# calling predictions without specifying new locations (newlocs) or 
# new covariates (newX) gives estimates of SVC only at the training locations.
pred.SVC <- predict(fit)

## ----visualization of prediction, fig.width=7, fig.height=7--------------
colnames(pred.SVC)[1:p] <- paste0("pred.",colnames(pred.SVC)[1:p])
coordinates(pred.SVC) <- ~loc_x+loc_y
all.SVC <- cbind(pred.SVC, sp.SVC[, 1:3])

# compute errors
all.SVC$err.SVC_1 <- all.SVC$pred.SVC_1 - all.SVC$SVC_1
all.SVC$err.SVC_2 <- all.SVC$pred.SVC_2 - all.SVC$SVC_2
all.SVC$err.SVC_3 <- all.SVC$pred.SVC_3 - all.SVC$SVC_3

colnames(all.SVC@data) <- paste0(rep(c("pred.", "true.", "err."), each = p), "SVC_", rep(1:p, 3))

spplot(all.SVC[, paste0(rep(c("true.", "err.", "pred."), each = p), 
                        "SVC_", 1:p)], colorkey = TRUE)

## ----n2500 figure, fig.width=7, fig.height=7-----------------------------
knitr::include_graphics("figures/SVCs_result_n2500_p3.png")

## ----n2500 example, eval=FALSE-------------------------------------------
#  # new m
#  m <- 50
#  
#  # new SVC model
#  sp.SVC <- varycoef:::fullSVC_reggrid(m = m, p = p,
#                                       pars = list(pars = pars,
#                                                   nugget.var = nugget.var))
#  spplot(sp.SVC, colorkey = TRUE)
#  
#  # total number of observations
#  n <- m^2
#  X <- matrix(c(rep(1, n), rnorm((p-1)*n)), ncol = p)
#  y <- apply(X * as.matrix(sp.SVC@data[, 1:p]), 1, sum) + sp.SVC@data[, p+1]
#  
#  
#  fit <- SVC_mle(y = y, X = X, locs = coordinates(sp.SVC))
#  
#  
#  
#  sp2500 <- predict(fit)
#  
#  
#  colnames(sp2500)[1:p] <- paste0("pred.",colnames(sp2500)[1:p])
#  coordinates(sp2500) <- ~loc_x+loc_y
#  all.SVC <- cbind(sp2500, sp.SVC[, 1:3])
#  
#  # compute errors
#  all.SVC$err.SVC_1 <- all.SVC$pred.SVC_1 - all.SVC$SVC_1
#  all.SVC$err.SVC_2 <- all.SVC$pred.SVC_2 - all.SVC$SVC_2
#  all.SVC$err.SVC_3 <- all.SVC$pred.SVC_3 - all.SVC$SVC_3
#  
#  colnames(all.SVC@data) <- paste0(rep(c("pred.", "true.", "err."), each = p), "SVC_", rep(1:p, 3))
#  
#  png(filename = "figures/SVCs_result_n2500_p3.png")
#  spplot(all.SVC[, paste0(rep(c("true.", "err.", "pred."), each = p),
#                          "SVC_", 1:p)], colorkey = TRUE,
#         as.table = TRUE, layout = c(3, 3))
#  dev.off()

## ----sample SVC model----------------------------------------------------
# new m
m <- 20

# number of fixed effects
pX <- 4
# number of mean-zero SVC
pW <- 3


# mean values
mu <- 1:pX

# new SVC model
sp.SVC <- varycoef:::fullSVC_reggrid(m = m, p = pW, 
                                     pars = list(pars = pars,
                                                 nugget.var = nugget.var), 
                                     seed = 4)

# total number of observations
n <- m^2
X <- matrix(c(rep(1, n), rnorm((pX-1)*n)), ncol = pX)
W <- X[, 1:pW]
# calculate y 
y <- 
  # X * mu
  X %*% mu +
  # W * beta_tilde
  apply(W * as.matrix(sp.SVC@data[, 1:pW]), 1, sum) + 
  # nugget
  sp.SVC@data[, pW+1]


# 
control <- SVC_mle.control()
control$init <- c(init[1:(2*pW + 1)], mu)
fit <- SVC_mle(y = y, X = X, W = W, locs = coordinates(sp.SVC), control = control)


## ----give predictions----------------------------------------------------
set.seed(5)
newlocs <- matrix(c(0, 0), ncol = 2)
X <- matrix(c(1, rnorm(pX-1)), ncol = pX)
W <- X[, 1:pW]

pred.SVC <- predict(fit, newlocs = newlocs)


ind.next <- apply((coordinates(sp.SVC) == 0) , 1, all)

sp.SVC[ind.next, ]

pred.SVC

