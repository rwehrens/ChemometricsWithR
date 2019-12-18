## Demo featuring all R code from Chapter 3 in the book "Chemometrics
## with R" by Ron Wehrens. The only additions to the book chapter code
## are loading the required libraries and data, and occasionally
## closing a graphics window.

graphics.off()
opar <- par(no.readonly = TRUE)

data(Prostate2000Raw)
prostate <- rowsum(t(Prostate2000Raw$intensity), 
                     group = rep(1:327, each = 2),
                     reorder = FALSE) / 2
dim(prostate)
prostate.type <- Prostate2000Raw$type[seq(1, 654, by = 2)]
matplot(Prostate2000Raw$mz[1:250],
        Prostate2000Raw$intensity[1:250, 1:2],
        type = "l", col = "pink", lty = 1,
        xlab = bquote(italic(.("m/z"))~.("(Da)")),
        ylab = "response")
lines(Prostate2000Raw$mz[1:250], prostate[1, 1:250], lwd = 2)

par(mfrow = c(1, 2))
rmeans <- rowMeans(embed(prostate[1, 1:250], 5))
plot(Prostate2000Raw$mz[1:250], prostate[1, 1:250],
     type = "l", xlab = "m/z", ylab = "response",
     main = "running means", col = 2)
lines(Prostate2000Raw$mz[3:248], rmeans, type = "l", lwd = 2)
plot(Prostate2000Raw$mz[1:250], prostate[1, 1:250],
     type = "l", xlab = "m/z", ylab = "response",
     main = "running median", col = 2)
lines(Prostate2000Raw$mz[1:250],
      runmed(prostate[1, 1:250], k = 5), type = "l", lwd = 2)

par(mfrow = c(1, 2))
plot(Prostate2000Raw$mz[-c(1, 2, 249:length(Prostate2000Raw$mz))],
     prostate[1, -c(1, 2, 249:length(Prostate2000Raw$mz))] -
     rmeans, ylim = c(-.2, .2),
     type = "l", xlab = "mz", ylab = "running mean: resids")
plot(Prostate2000Raw$mz[1:250],
     prostate[1, 1:250] - runmed(prostate[1, 1:250], k = 5),
     ylim = c(-.2, .2),
     type = "l", xlab = "mz", ylab = "running median: resids")

par(opar)
mznew <- colMeans(matrix(Prostate2000Raw$mz[1:250], nrow = 5))
xnew <- colMeans(matrix(prostate[1, 1:250], nrow = 5))
plot(Prostate2000Raw$mz[1:250], prostate[1, 1:250], 
     type = "l", xlab = "m/z", ylab = "response",
     main = "binning", col = 2)
lines(mznew, xnew, lwd = 2)

## ########################################################################
data(gasoline, package = "pls")
wavelengths <- seq(900, 1700, by = 2)

require(signal)
nir.deriv <- apply(gasoline$NIR, 1, sgolayfilt, m = 1)
nir.diff <- t(apply(gasoline$NIR, 1, diff))
matplot(wavelengths[-1] + 1, t(nir.diff),
        xlab = "Wavelength (nm)", ylab = "1/R (1st deriv.)",
        type = "n")
abline(h = 0, col = "gray")
matlines(wavelengths[-1] + 1, t(nir.diff), lty = 1)

require(pls)
nir.msc <- msc(gasoline$NIR)

require(ptw)
par(mfrow = c(1, 2))
x <- Prostate2000Raw$intensity[1:4000, 1]
mz <- Prostate2000Raw$mz[1:4000]
lsection <- 200
xmat <- matrix(x, nrow = lsection)
ymin <- apply(xmat, 2, min)
plot(mz, x, type = "l", col = "darkgray", ylim = c(-1.2, 5),
     xlab = "m/z", ylab = "I")
lines(mz, rep(ymin, each = lsection))
bsln.loess <- loess(ymin ~ mz[seq(101, 4000, by = 200)])
lines(mz, predict(bsln.loess, mz), lwd = 2, col = "blue")
plot(mz, x, col = "darkgray", type = "l", ylim = c(-1.2, 5),
     xlab = "m/z", ylab = "I")
lines(mz, asysm(x), lwd = 2, col = "blue")
par(opar)

data(lcms)
plot(time, lcms[1, , 2], type = "l",
     xlab = "Time (s)", ylab = "I", main = "Mass channel 1")
lines(time, lcms[1, , 3], type = "l", col = 2)
legend("topleft", legend = paste("Sample", 2:3),
       lty = 1, col = 1:2)

sref <- lcms[1, , 2]
ssamp <- lcms[1, , 3]
lcms.warp <- ptw(sref, ssamp, init.coef = c(0, 1, 0))
summary(lcms.warp)

plot(time, sref, type = "l", lwd = 2, col = 2,
     xlim = time[c(600, 1300)], xlab = "Time (s)", ylab = "I")
lines(time, ssamp + 1e6, lty = 2)
lines(time, lcms.warp$warped.sample + 2e6, col = 4)
legend("topleft", lty = c(1, 2, 1), col = c(2, 1, 4),
       legend = c("Reference", "Sample", "Warped sample"),
       lwd = c(2, 1, 1))

lcms.warpRMS <- ptw(sref, ssamp, optim.crit = "RMS")
lcms.warpRMS$warp.coef

lcms.warp2 <- ptw(sref, ssamp, init = c(0, 1, 0, 0))
lcms.warp3 <- ptw(sref, ssamp, init = c(0, 1, 0, 0, 0))
lcms.warp4 <- ptw(sref, ssamp, init = c(0, 1, 0, 0, 0, 0))

allwarps <- list(lcms.warp, lcms.warp2, lcms.warp3, lcms.warp4)
wccs <- sapply(allwarps, function(x) x$crit.value)

allwarp.funs <- sapply(allwarps, function(x) x$warp.fun)
warpings <- allwarp.funs - 1:length(sref)

matplot(time, warpings, type = "l", lty = rep(c(1, 2), 2), 
        col = 1:4, ylim = c(min(warpings), 0.5))
abline(h = 0, col = "gray", lty = 2)
legend("bottom", lty = rep(c(1, 2), 2), col = 1:4, bty = "n",
       legend = paste("Degree", 2:5, " - WCC =", round(wccs, 3)))

srefM <- lcms[, , 2]
ssampM <- lcms[, , 3]
traces <- select.traces(srefM, criterion = "var")
lcms.warpglobal <-
  ptw(srefM, ssampM, warp.type = "global",
      selected.traces = traces$trace.nrs[1:10])

summary(lcms.warpglobal)

npad <- 500
srefM2 <- padzeros(lcms[, , 2], npad, "both")
ssampM2 <- padzeros(lcms[, , 3], npad, "both")
sample2.indiv.warp <-
  ptw(srefM2, ssampM2, init.coef = c(0, 1, 0, 0, 0, 0, 0, 0, 0))
sample2.global.warp <-
  ptw(srefM2, ssampM2, init.coef = c(0, 1, 0),
      warp.type = "global")
sample2.final.warp <-
  ptw(srefM2, ssampM2,
      init.coef = c(sample2.global.warp$warp.coef, 0, 0))

par(mfrow = c(2, 1))
orig.ind <- (npad + 1):(npad + 2000)
plot(time, colSums(srefM2)[orig.ind], col = 2, 
     type = "l", main = "PTW (indiv.)", ylab = "I")
lines(time, colSums(ssampM2)[orig.ind], col = "gray")
lines(time,
      colSums(sample2.indiv.warp$warped.sample)[orig.ind], lwd = 2)
legend("topleft", bty = "n", lty = 1, lwd = c(1, 2, 1),
       legend = c("Original sample", "Warped sample", "Reference"),
       col = c("gray", "black", "red"))
plot(time, colSums(srefM2)[orig.ind], col = 2, 
     type = "l", main = "PTW (global + indiv.)", ylab = "I")
lines(time, colSums(ssampM2)[orig.ind], col = "gray")
lines(time, colSums(sample2.final.warp$warped.sample)[orig.ind], lwd = 2)
legend("topleft", bty = "n", lty = 1, lwd = c(1, 2, 1),
       legend = c("Original sample", "Warped sample", "Reference"),
       col = c("gray", "black", "red"))

require(dtw, quiet = TRUE)
sample2.dtw <- matrix(0, 100, 2000)
for (i in 1:100) {                 
  warpfun.dtw.i <- dtw(ssampM[i, ], srefM[i, ], keep = TRUE)
  new.indices <- warp(warpfun.dtw.i, index.reference = FALSE)
  sample2.dtw[i, ] <- ssampM[i, new.indices]
}                                                                        

par(mfrow = c(1, 2))
warpfun.dtw <- dtw(ssamp, sref)
plot(warpfun.dtw)
abline(0, 1, col = "gray", lty = 2)
par(pty = "s")
warpfun.dtw <- dtw(ssamp, sref, keep.internals = TRUE)
contour(warpfun.dtw$costMatrix,
        x = 1:2000, y = 1:2000,
        drawlabels = FALSE,
        xlab = "Query index", ylab = "Reference index")
lines(warpfun.dtw$index1, warpfun.dtw$index2, col = 2, lwd = 3)

par(opar)
warped.signal <- warp(warpfun.dtw, index.reference = FALSE)
plot(time, sref, type = "l", lwd = 2, col = 2,
     xlim = c(time[600], time[1300]),
     xlab = "Time", ylab = "I")
lines(time, ssamp + 1e6, lty = 2)
lines(time, ssamp[warped.signal] + 2e6, col = 4)
legend("topleft", lty = c(1, 2, 1), col = c(2, 1, 4),
       legend = c("Reference", "Sample", "Warped sample"),
       lwd = c(2, 1, 1))

warp.dtw.gl <- dtw(t(ssampM), t(srefM))
samp.aligned <- lcms[, warp(warp.dtw.gl), 3]    

par(mfrow = c(2, 1))
plot(time, colSums(srefM), col = 2,
     type = "l", main = "DTW (indiv.)", ylab = "I")
lines(time, colSums(ssampM), col = "gray")
lines(time, colSums(sample2.dtw), lwd = 2)
legend("topleft", legend = c("Sample", "Reference"),
       lty = 1, col = 1:2, lwd = 2:1)

plot(time, colSums(srefM), col = 2,
     type = "l", main = "DTW (global)", ylab = "I")
lines(time, colSums(ssampM), col = "gray")
lines(time, colSums(samp.aligned), lwd = 2)
legend("topleft", legend = c("Sample", "Reference"),
       lty = 1, col = 1:2, lwd = 2:1)

pick.peaks <- function(x, span) {
  span.width <- span * 2 + 1
  loc.max <- span.width + 1 -
    apply(embed(x, span.width), 1, which.max)
  loc.max[loc.max == 1 | loc.max == span.width] <- NA

  pks <- loc.max + 0:(length(loc.max) - 1)
  unique(pks[!is.na(pks)])
}

par(mfrow = c(1, 2))
pks10 <- pick.peaks(rmeans, 10)
plot(Prostate2000Raw$mz[3:248], rmeans, type = "l",
     xlab = "m/z", ylab = "Response", main = "span = 10")
abline(v = Prostate2000Raw$mz[pks10 + 2], col = 2)
pks40 <- pick.peaks(rmeans, 40)
plot(Prostate2000Raw$mz[3:248], rmeans, type = "l",
     xlab = "m/z", ylab = "Response", main = "span = 40")
abline(v = Prostate2000Raw$mz[pks40 + 2], col = 2)

range(apply(prostate[1:10, ], 1, max))
range(rowSums(prostate[1:10, ]))

prost10.rangesc <- sweep(prostate[1:10, ], MARGIN = 1,
                         apply(prostate[1:10, ], 1, max),
                         FUN = "/")
apply(prost10.rangesc, 1, max)
range(rowSums(prost10.rangesc))

prost10.lengthsc <- sweep(prostate[1:10, ], MARGIN = 1,
                          apply(prostate[1:10, ], 1,
                                function(x) sqrt(sum(x^2))),
                          FUN = "/")
range(apply(prost10.lengthsc, 1, max))
range(rowSums(prost10.lengthsc))

prost10.varsc <- sweep(prostate[1:10, ], MARGIN = 1,
                       apply(prostate[1:10, ], 1, sd),
                       FUN = "/")
range(apply(prost10.varsc, 1, max))
range(rowSums(prost10.varsc))

NIR.mc <- t(sweep(gasoline$NIR, 2, colMeans(gasoline$NIR)))

par(opar)
NIR.mc <- scale(gasoline$NIR, scale = FALSE)
matplot(wavelengths, t(NIR.mc),
        type = "l", xlab = "Wavelength (nm)",
        ylab = "1/R (mean-centered)", lty = 1)

options("width" = 60)
data(wines, package = "kohonen")
apply(wines, 2, range)

par(mfrow = c(1, 2))
wines.mc <- scale(wines, scale = FALSE)
wines.sc <- scale(wines, scale = TRUE)
boxplot(wines.mc ~ col(wines.mc), 
        main = "Mean-centered wine data")
boxplot(wines.sc ~ col(wines.sc), 
        main = "Autoscaled wine data")

