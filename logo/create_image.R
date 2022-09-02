library(varycoef)

par(mar = rep(1, 4))
matplot(
  # Locations and SVCs
  x = SVCdata$locs, 
  y = SVCdata$beta, 
  # Lines
  type = "l", 
  lty = 1,
  lwd = 8, 
  # no labels (will follow), no box
  xlab = "", 
  ylab = "",
  bty = "n", 
  xaxt = "n", 
  yaxt = "n",
  col = 2:3
)
# Adding axis
axis(1, at = 2*(0:5), lwd=4, lwd.tick=4, lab=F, cex.axis = 1.5, las = 1)
axis(2, at = 0:3, lwd=4, lwd.tick=4, lab=F, cex.axis = 1.5, las = 1)
