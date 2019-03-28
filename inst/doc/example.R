## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(varycoef)

## ----define SVC----------------------------------------------------------
# number of SVC
p <- 3

(pars <- data.frame(mu = rep(0, p), 
                    var = c(0.1, 0.2, 0.3), 
                    scale = c(0.3, 0.1, 0.2)))
nugget.var <- 0.05

## ----sample SVC, fig.width=7, fig.height=7-------------------------------
library(RandomFields)
library(sp)
m <- 20

# number of observations
n <- m^2

# regular grid locations
locs <- expand.grid(x = seq(0, 1, length.out = m), 
                    y = seq(0, 1, length.out = m))

set.seed(123)

# SVC model
model <- apply(pars, 1, function(x) {
  RFsimulate(RMexp(x["var"], x["scale"]), 
             x = locs[, "x"], y = locs[, "y"])
})

model[[p+1]] <- RFsimulate(RMnugget(var = nugget.var), 
                           x = locs[, "x"], y = locs[, "y"])
sp.SVC <- Reduce(cbind, model)
sp.SVC <- SpatialPointsDataFrame(coords = sp.SVC@coords, 
                                 data = sp.SVC@data)
colnames(sp.SVC@data) <- c(paste0("SVC_", 1:p), "nugget")

spplot(sp.SVC, colorkey = TRUE)

## ----sample covariates---------------------------------------------------
X <- matrix(c(rep(1, n), rnorm((p-1)*n)), ncol = p)
head(X)

## ----compute response----------------------------------------------------
y <- apply(X * as.matrix(sp.SVC@data[, 1:p]), 1, sum) + sp.SVC@data[, p+1]

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



lower <- c(rep(0, 2*p+1), rep(-Inf, p))

## ----SVC MLE-------------------------------------------------------------
fit <- SVC_mle(y = y, X = X, locs = locs, init = init, 
               lower = lower)

class(fit)

# comparison of estimated and true parameters
rbind(fit$optim.output$par, 
      c(pars[, "scale"], pars[, "var"], nugget.var, pars[, "mu"]))

## ----make predictions----------------------------------------------------
# calling predictions without specifying new locations (newlocs) or 
# new covariates (newX) gives estimates of SVC only at the training location.
pred.SVC <- predict(fit)

## ----visualization of prediction, fig.width=7, fig.height=7--------------
colnames(pred.SVC)[1:p] <- paste0("pred.",colnames(pred.SVC)[1:p])
coordinates(pred.SVC) <- ~loc_x+loc_y
all.SVC <- cbind(pred.SVC, sp.SVC[, 1:3])

# compute errors
all.SVC$err.SVC_1 <- all.SVC$pred.SVC_1 - all.SVC$SVC_1
all.SVC$err.SVC_2 <- all.SVC$pred.SVC_2 - all.SVC$SVC_2
all.SVC$err.SVC_3 <- all.SVC$pred.SVC_3 - all.SVC$SVC_3


spplot(all.SVC[, paste0(rep(c("pred.", "err."), each = p), 
                        "SVC_", 1:p)], colorkey = TRUE)

## ----n1600 example, fig.width=7, fig.height=7, echo = FALSE--------------
knitr::include_graphics("figures/SVCs_result_n1600_p3.png")

