## Demo featuring all R code from Chapter 5 in the book "Chemometrics
## with R" by Ron Wehrens. The only additions to the book chapter code
## are loading the required libraries and data, and occasionally
## closing a graphics window.
graphics.off()
opar <- par(no.readonly = TRUE)

require(kohonen, quiet = TRUE)

require(MASS, quiet = TRUE) # for mvrnorm
set.seed(17)
rpoints <- rbind(mvrnorm(100, rep(0, 2), diag(2)),
                 mvrnorm(500, rep(0, 2), diag(2) * .025))
rpoints.som <- som(rpoints, somgrid(5, 5, "rectangular"))
par(mfrow = c(1, 2))
plot(rpoints, col = "blue", xlab = "", ylab = "", axes = FALSE)
X <- getCodes(rpoints.som, 1)
points(X)
text(X[c(1, 5, 21, 25), ], labels = c(1, 5, 21, 25), pos = c(1, 3, 1, 3))
horiz1 <- (1:25)[-seq(1, 25, by = 5)]
vertic1 <- 6:25
segments(X[horiz1, 1], X[horiz1, 2],
         X[horiz1 - 1, 1], X[horiz1 - 1, 2])
segments(X[vertic1, 1], X[vertic1, 2],
         X[vertic1 - 5, 1], X[vertic1 - 5, 2])

plot(rpoints.som, type = "mapping", col = "blue", pch = 1,
     main = "", border = "transparent")
par("col" = "black")
X <- rpoints.som$grid$pts
points(X)
text(X[c(1, 5, 21, 25), ], labels = c(1, 5, 21, 25), pos = c(1, 1, 3, 3))
horiz1 <- (1:25)[-seq(1, 25, by = 5)]
vertic1 <- 6:25
segments(X[horiz1, 1], X[horiz1, 2],                                           
         X[horiz1 - 1, 1], X[horiz1 - 1, 2])                                   
segments(X[vertic1, 1], X[vertic1, 2],                                         
         X[vertic1 - 5, 1], X[vertic1 - 5, 2])

data(wines)
wines.sc <- scale(wines)
set.seed(7)
wines.som <- som(wines.sc, somgrid(5, 4, "hexagonal"))
par(opar)
plot(wines.som, type = "codes")

par(mfrow = c(2, 2))
for (i in c(1, 8, 11, 13)) 
  plot(wines.som, "property",
       property = getCodes(wines.som, 1)[, i],
       main = colnames(wines)[i])

par(opar)
plot(wines.som, type = "mapping",
     col = as.integer(vintages), pch = as.integer(vintages))

plot(wines.som, type = "dist.neighb")

summary(wines.som)

par(mfrow = c(1, 2))
plot(wines.som, "changes")
plot(wines.som, "quality")

data(Prostate2000Raw)
set.seed(7)
X <- t(Prostate2000Raw$intensity)
prostate.som <- som(X, somgrid(7, 5, "hexagonal"))

par(opar)
require(lattice, quiet = TRUE)
types <- as.integer(Prostate2000Raw$type)
trellis.cols <- trellis.par.get("superpose.symbol")$col[c(2, 3, 1)]
plot(prostate.som, "mapping", col = trellis.cols[types],
     pch =  types, main = "Prostate data")
legend("bottom", legend = levels(Prostate2000Raw$type),
       col = trellis.cols, pch = 1:3, ncol = 3, bty = "n")

units <- c(7, 21, 35)
unitfs <- paste("Unit", units)
prost.plotdf <-
  data.frame(mz = Prostate2000Raw$mz,
             intensity = c(t(getCodes(prostate.som, 1)[units, ])),
             unit = rep(factor(unitfs, levels = unitfs),
                        each = length(Prostate2000Raw$mz)))
xyplot(intensity ~ mz | unit, data = prost.plotdf, type = "l",
       scale = list(y = "free"), as.table = TRUE,
       xlab = bquote(italic(.("m/z"))~.("(Da)")),
       groups = unit, layout = c(1, 3),
       panel = function(...) {
         panel.abline(v = c(3300, 4000, 6000, 6200),
                      col = "gray", lty = 2)
         panel.xyplot(...)
       })

