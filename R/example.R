fullSVC_reggrid <- function(m = 10, p = 3, pars = NULL, seed = 123) {

  # number of observations
  n <- m^2

  # regular grid locations
  locs <- expand.grid(x = seq(0, 1, length.out = m),
                      y = seq(0, 1, length.out = m))

  set.seed(seed)

  # SVC model
  model <- apply(pars$pars, 1, function(x) {
    RandomFields::RFsimulate(
      RandomFields::RMexp(x["var"], x["scale"]),
                          x = locs[, "x"], y = locs[, "y"])
  })

  model[[p+1]] <- RandomFields::RFsimulate(
    RandomFields::RMnugget(var = pars$nugget.var),
                           x = locs[, "x"], y = locs[, "y"])
  sp.SVC <- Reduce(cbind, model)
  sp.SVC <- sp::SpatialPointsDataFrame(coords = sp.SVC@coords,
                                       data = sp.SVC@data)
  colnames(sp.SVC@data) <- c(paste0("SVC_", 1:p), "nugget")

  return(sp.SVC)
}
