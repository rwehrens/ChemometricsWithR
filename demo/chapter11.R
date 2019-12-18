## Demo featuring all R code from Chapter 11 in the book "Chemometrics
## with R" by Ron Wehrens. The only additions to the book chapter code
## are loading the required libraries and data, and occasionally
## closing a graphics window.

graphics.off()
opar <- par(no.readonly = TRUE)

data(arabidopsis)
naLimitPerc <- 40
naLimit <- floor(nrow(arabidopsis) * naLimitPerc / 100)
nNA <- apply(arabidopsis, 2, function(x) sum(is.na(x)))
naIdx <- which(nNA < naLimit)
X.ara <- arabidopsis[, naIdx]

plot(sort(nNA), type = "h", col = 4, xlab = "Variables (sorted)",
     ylab = "Number of missing values per metabolite",
     main = "Arabidopsis data")
abline(h = naLimit, col = "darkgray")
text(1, naLimit, naLimit, pos = 3)
abline(v = length(naIdx) + .5, col = "darkgray", lty = 2)
text(length(naIdx), max(nNA), length(naIdx), pos = 4)

X.ara.l <- log(X.ara)
X.ara.l.sc <- scale(X.ara.l)

X.ara.cov <- cov(t(X.ara.l.sc), use = "pairwise.complete.obs")
sum(is.na(X.ara.cov))
X.ara.svd <- svd(X.ara.cov, nu = 2, nv = 2)
ara.PCA.svd <-
  structure(
    list(scores = X.ara.svd$u %*% diag(sqrt(X.ara.svd$d[1:2])),
         var = X.ara.svd$d,
         totalvar = sum(X.ara.svd$d),
         centered.data = TRUE),
    class = "PCA")

X.ara.imput1 <-
  apply(X.ara, 2,
        function(x) {
          x[is.na(x)] <- min(x, na.rm = TRUE)
          x
          })

require(missMDA, quiet = TRUE)
X.ara.pcaimput <- imputePCA(X.ara.l, ncp = 2)$completeObs

par(mfrow = c(1, 3))
scoreplot(ara.PCA.svd, main = "PCA using cov")
ara.PCA.minimputation <- PCA(scale(log(X.ara.imput1)))
scoreplot(ara.PCA.minimputation,
          main = "PCA using minimum imputation")
ara.PCA.pcaimputation <- PCA(scale(X.ara.pcaimput))
scoreplot(ara.PCA.pcaimputation,
          main = "PCA using PCA imputation")

X.ara.imput1 <-
  apply(X.ara, 2,
        function(x) {
          x[is.na(x)] <- min(x, na.rm = TRUE)
          x
          })

par(opar)
p1 <- hist(X.ara.l, plot = FALSE)
p2 <- hist(X.ara.pcaimput[is.na(X.ara.l)],
           plot = FALSE)
p3 <- hist(log(X.ara.imput1)[is.na(X.ara.l)], plot = FALSE,
           breaks = p2$breaks)
plot(p1, xlab = "Feature intensity (log scale)",
     main = "Arabidopsis data")
col3 <- rgb(.5, .5, 0, .2)
col2 <- rgb(0, .5, .5, .2)
plot(p3, col = col3, add = TRUE)
plot(p2, col = col2, add = TRUE)

legend("topright",
       legend = c("Measured values", "Imputed values (min)",
                  "Imputed values (PCA)"),
       fill = c("white", col3, col2), bty = "n")

rownames(X.ara.l) <- rep("", nrow(X.ara.l))
colnames(X.ara.l) <- paste("V", 1:ncol(X.ara.l), sep = "")
countNAs <- apply(X.ara.l[, 1:20], 2, function(x) sum(is.na(x)))
countNAs[countNAs > 0]

ara.PCA.Minput <- MIPCA(X.ara.l[, 1:20], ncp = 2, scale = TRUE)

par(mfrow = c(1, 2))
plot(ara.PCA.Minput, choice = "dim", new.plot = FALSE)
plot(ara.PCA.Minput, choice = "var", new.plot = FALSE)

detach("package:missMDA", unload = TRUE)

require(MASS, quiet = TRUE)

data(wines, package = "kohonen")
wine.classes <- as.integer(vintages)
wines.sc <- scale(wines)
wines.odd <- seq(1, nrow(wines), by = 2)
wines.even <- seq(2, nrow(wines), by = 2)
twowines <- wines[vintages != "Barolo",]
twovintages <- factor(vintages[vintages != "Barolo"])

par(mfrow = c(1, 2))
X <- wines[vintages == "Grignolino", ]
X.sc <- scale(X)
X.clPCA <- princomp(X.sc)
X.robPCA <- princomp(X.sc, covmat = cov.mcd(X.sc))
biplot(X.clPCA, main = "Classical PCA")
biplot(X.robPCA, main = "MCD-based PCA")

require(rrcov, quiet = TRUE)
X.HubPCA5 <- PcaHubert(X.sc, k = 5)
summary(X.HubPCA5)

X.HubPCA <- PcaHubert(X.sc)
summary(X.HubPCA)

par(opar)
plot(X.HubPCA)

require(qcc, quiet = TRUE)
metabNr <- 2
B1 <- which(arabidopsis.Y$Batch == "B1")
B23 <- which(arabidopsis.Y$Batch %in% c("B2", "B3"))
ara.qcc <- qcc(data = arabidopsis[B1, metabNr], type = "xbar.one",
               newdata = arabidopsis[B23, metabNr],
               plot = FALSE)
ara.cusum <- cusum(data = arabidopsis[B1, metabNr],
                   newdata = arabidopsis[B23, metabNr],
                   plot = FALSE)

par(mfrow = c(1, 2))
qcc.options(bg.margin = "transparent",
            beyond.limits = list(pch = 15, col = "red"),
            violating.runs = list(pch = 17, col = "orange"))
ara.qcc$data.name <- ara.cusum$data.name <- "Batch 1"
ara.qcc$newdata.name <- ara.cusum$newdata.name <- "Batches 2, 3"
names(ara.qcc$statistics) <- names(ara.cusum$statistics) <-
  names(ara.qcc$newstats) <- names(ara.cusum$newstats) <- " "
plot(ara.qcc, add.stats = FALSE, restore.par = FALSE)
plot(ara.cusum, add.stats = FALSE)

require(MSQC, quiet = TRUE)
par(mfrow = c(1, 2))
idx <- which(apply(arabidopsis, 2, function(x) sum(is.na(x)) == 0))
chi2chart <-
  capture.output(mult.chart(arabidopsis[c(B1, B23), idx], type = "chi",
                            Xmv = colMeans(arabidopsis[B1, idx]),
                            S = cov(arabidopsis[B1, idx])))
abline(v = length(B1) + .5, lty = 2)

MCUSUMchart <- mult.chart(arabidopsis[c(B1, B23), idx], type = "mcusum2",
                          Xmv = colMeans(arabidopsis[B1, idx]),
                          S = cov(arabidopsis[B1, idx]))
abline(v = length(B1) + .5, lty = 2)

require(pls, quiet = TRUE)
data(gasoline, package = "pls")
wavelengths <- seq(900, 1700, by = 2)
gas.odd <- seq(1, nrow(gasoline$NIR), by = 2)
gas.even <- seq(2, nrow(gasoline$NIR), by = 2)

gasolineSC <- gasoline
gasolineSC$NIR <-
  scale(gasolineSC$NIR, scale = FALSE,
        center = colMeans(gasolineSC$NIR[gas.odd, ]))
gasolineSC.pls <- plsr(octane ~ ., data = gasolineSC, ncomp = 5,
                       subset = gas.odd, validation = "LOO")
ww <- gasolineSC.pls$loading.weights[, 1]
pp <- gasolineSC.pls$loadings[, 1]
w.ortho <- pp - c(crossprod(ww, pp)/crossprod(ww)) * ww
t.ortho <- gasolineSC$NIR[gas.odd, ] %*% w.ortho
p.ortho <- crossprod(gasolineSC$NIR[gas.odd, ], t.ortho) /
  c(crossprod(t.ortho))
Xcorr <- gasolineSC$NIR[gas.odd, ] - tcrossprod(t.ortho, p.ortho)

gasolineSC.osc1 <- data.frame(octane = gasolineSC$octane[gas.odd],
                              NIR = Xcorr)
gasolineSC.opls1 <- plsr(octane ~ ., data = gasolineSC.osc1,
                         ncomp = 5, validation = "LOO")

pp2 <- gasolineSC.opls1$loadings[, 1]
w.ortho2 <- pp2 - c(crossprod(ww, pp2)/crossprod(ww)) * ww
t.ortho2 <- Xcorr %*% w.ortho2
p.ortho2 <- crossprod(Xcorr, t.ortho2) / c(crossprod(t.ortho2))
Xcorr2 <- Xcorr - tcrossprod(t.ortho2, p.ortho2)
gasolineSC.osc2 <- data.frame(octane = gasolineSC$octane[gas.odd],
                            NIR = Xcorr2)
gasolineSC.opls2 <- plsr(octane ~ ., data = gasolineSC.osc2,
                       ncomp = 5, validation = "LOO")

par(opar)
plot(gasolineSC.pls, "validation", estimate = "CV",
     ylim = c(0.2, 1.5),
     main = "GasolineSC training data (validation)")
lines(0:5, c(RMSEP(gasolineSC.opls1, estimate = "CV"))$val,
      col = 2, lty = 2)
lines(0:5, c(RMSEP(gasolineSC.opls2, estimate = "CV"))$val,
      col = 4, lty = 4)
legend("topright", lty = c(1, 2, 4), col = c(1, 2, 4),
       legend = c("PLS", "OPLS: 1 OSC component", 
                  "OPLS: 2 OSC components"))

Xtst <- gasolineSC$NIR[gas.even, ]
t.tst <- Xtst %*% w.ortho
p.tst <- crossprod(Xtst, t.tst) / c(crossprod(t.tst))
Xtst.osc1 <- Xtst - tcrossprod(t.tst, p.tst)
gasolineSC.opls1.pred <- predict(gasolineSC.opls1,
                               newdata = Xtst.osc1,
                               ncomp = 2)

t.tst2 <- Xtst.osc1 %*% w.ortho2
p.tst2 <- crossprod(Xtst.osc1, t.tst2) / c(crossprod(t.tst2))
Xtst.osc2 <- Xtst.osc1 - tcrossprod(t.tst2, p.tst2)
gasolineSC.opls2.pred <- predict(gasolineSC.opls2,
                               newdata = Xtst.osc2,
                               ncomp = 1)

RMSEP(gasolineSC.pls, newdata = gasolineSC[gas.even, ],
      ncomp = 3, intercept = FALSE)
rms(gasolineSC$octane[gas.even], gasolineSC.opls1.pred)
rms(gasolineSC$octane[gas.even], gasolineSC.opls2.pred)

require(glmnet, quiet = TRUE)
set.seed(7)
data(spikedApples, package = "BioMark")
X <- sqrt(spikedApples$dataMatrix)
Y <- rep(0:1, each = 10)
apple.df <- data.frame(Y = Y, X = X)
apple.pls <- plsr(Y ~ X, data = apple.df, ncomp = 5,
                  validation = "LOO")
nPLS <- selectNcomp(apple.pls, method = "onesigma")
apple.lasso <- cv.glmnet(X, Y, family = "binomial")

tvals <- apply(X, 2,
               function(x) t.test(x[1:10], x[11:20])$statistic)
allcoefs <-
  data.frame(studentt = tvals,
             pls = c(coef(apple.pls, ncomp = nPLS)),
             lasso = coef(apple.lasso,
                          lambda = apple.lasso$lambda.1se)[-1])

par(opar)
N <- ncol(spikedApples$dataMatrix)
biom <- spikedApples$biom
nobiom <- (1:N)[-spikedApples$biom]
pairs(allcoefs,
      panel = function(x, y, ...) {
        abline(h = 0, v = 0, col = "lightgray", lty = 2)
        points(x[nobiom], y[nobiom], col = "lightgray")
        points(x[biom], y[biom], cex = 2)
      })

require(BioMark, quiet = TRUE)
biomarkerSets <-
  get.biom(X, factor(Y), type = "coef", ncomp = nPLS,
           fmethod = c("studentt", "pls", "lasso"))

set.seed(17)
apple.stab <- get.biom(X = X, Y = factor(Y), ncomp = 1:3,
                       type = "stab", fmethod = c("pls", "pcr"))
selected.variables <- selection(apple.stab)
unlist(sapply(selected.variables, function(x) sapply(x, length)))

data(shootout)
wl <- seq(600, 1898, by = 2)
indices <- which(wl >= 1100 & wl <= 1700)
nir.training1 <-
  data.frame(X = I(shootout$calibrate.1[, indices]), 
             y = shootout$calibrate.Y[, 3])
mod1 <- plsr(y ~ X, data = nir.training1,
             ncomp = 5, validation = "LOO")
RMSEP(mod1, estimate = "CV")

nir.training2 <-
  data.frame(X = I(shootout$calibrate.2[, indices]), 
             y = shootout$calibrate.Y[, 3])
mod2 <- plsr(y ~ X, data = nir.training2,
ncomp = 5, validation = "LOO")

par(opar)
layout(matrix(c(1, 1, 1, 2, 2), nrow = 1))             
plot(seq(1100, 1700, by = 2), coef(mod1, ncomp = 3), type = "l",
     xlab = "wavelength (nm)", ylab = "model coefficients",
     col = 1)
lines(seq(1100, 1700, by = 2), coef(mod2, ncomp = 3), col = "red",
      lty = 2, lwd = 2)
legend("top", legend = c("set 1", "set 2"), bty = "n",
       col = 1:2, lty = 1:2, lwd = 1:2)
plot(seq(1100, 1700, by = 2), coef(mod1, ncomp = 3), type = "l",
     xlab = "wavelength (nm)", ylab = "model coefficients",
     col = 1, xlim = c(1600, 1700))
lines(seq(1100, 1700, by = 2), coef(mod2, ncomp = 3), col = "red",
      lty = 2, lwd = 2)

RMSEP(mod1, estimate = "test", ncomp = 3, intercept = FALSE,
      newdata = data.frame(y = shootout$test.Y[, 3],
                           X = I(shootout$test.1[, indices])))

RMSEP(mod1, estimate = "test", ncomp = 3, intercept = FALSE,
      newdata = data.frame(y = shootout$test.Y[, 3],
                           X = I(shootout$test.2[, indices])))

recalib.indices <- 1:5 * 10
F1 <- ginv(shootout$calibrate.2[recalib.indices, indices]) %*%
  shootout$calibrate.1[recalib.indices, indices]
RMSEP(mod1, estimate = "test", ncomp = 3, intercept = FALSE,
      newdata = data.frame(y = shootout$test.Y[, 3], 
        X = I(shootout$test.2[, indices] %*% F1)))

require(lattice, quiet = TRUE)
levelplot(F1, contour = TRUE)

require(fpc, quiet = TRUE)
aveBhatta <- function(Xmat, batches) {
  N <- nlevels(batches)
  batch.means <- lapply(levels(batches),
                        function(btch)
                          colMeans(Xmat[batches == btch, ]))
  batch.covs <- lapply(levels(batches),
                       function(btch)
                         cov(Xmat[batches == btch, ]))

  alldists <- 0
  for (ii in 1:(N-1))
    for (jj in (ii+1):N)
      alldists <- alldists +
            bhattacharyya.dist(batch.means[[jj]], batch.means[[ii]],
                               batch.covs[[jj]], batch.covs[[ii]])

  alldists / (0.5*N*(N-1))
}

par(opar)
Xsample <- arabidopsis
Xsample <- Xsample[, apply(Xsample, 2, function(x) !all(is.na(x)))]
for (i in 1:ncol(Xsample))
  Xsample[is.na(Xsample[, i]), i] <- mean(Xsample[, i], na.rm = TRUE)
Xsample <- Xsample[, apply(Xsample, 2, sd, na.rm = TRUE) > 0]
X.PCA <- PCA(scale(Xsample))
pchs  <- as.integer(arabidopsis.Y$Batch)
cols <- rep("gray", nrow(Xsample))
cols[arabidopsis.Y$Batch  == "B1"] <- "red"
cols[arabidopsis.Y$Batch  == "B2"] <- "blue"
cols[arabidopsis.Y$Batch  == "B8"] <- "purple"
lcols <- rep("gray", 10)
lcols[c(1, 2, 8)] <- c("red", "blue", "purple")
scoreplot(X.PCA, col = cols, pch = pchs)
legend("bottomright", legend = levels(arabidopsis.Y$Batch),
       ncol = 2, pch = 1:nlevels(arabidopsis.Y$Batch), col = lcols)

bhatta.orig <- aveBhatta(scores(X.PCA, 2), arabidopsis.Y$Batch)
title(main = paste("Average between-batch distance:",
                   round(bhatta.orig, 3)))

ara.df <- cbind(data.frame(X = arabidopsis[, 2]),
                arabidopsis.Y)
ref.idx <- ara.df$SCode == "ref"
plot(X ~ SeqNr, data = ara.df, ylim = c(20, 24),
     col = as.numeric(ref.idx) + 1,
     pch = c(1, 19)[as.numeric(ref.idx) + 1],
     xlab = "Injection number", ylab = "Intensity (log-scaled)",
     main = paste("Metabolite 2 before correction"))
batch.lims <- aggregate(ara.df$SeqNr,
                        by = list(ara.df$Batch),
                        FUN = range)$x
abline(v = batch.lims[-1, 1] - 0.5, lty = 2)
BLines <- lm(X ~ SeqNr * Batch, data = ara.df,
             subset = SCode == "ref")
ara.df$QCpredictions <- predict(BLines, newdata = ara.df)
for (ii in levels(ara.df$Batch))
  lines(ara.df$SeqNr[ara.df$Batch == ii],
        predict(BLines, newdata = ara.df[ara.df$Batch == ii, ]),
        col = 2)
ara.df$corrected <-
  ara.df$X - ara.df$QCpredictions + mean(ara.df$X)
plot(corrected ~ SeqNr, data = ara.df, ylim = c(20, 24),
     col = as.numeric(ref.idx) + 1,
     pch = c(1, 19)[as.numeric(ref.idx) + 1],
     xlab = "Injection number", ylab = "Intensity (log-scaled)",
     main = paste("Metabolite 2 after correction"))
abline(v = batch.lims[-1, 1] - 0.5, lty = 2)
BLines2 <- lm(corrected ~ SeqNr * Batch, data = ara.df,
              subset = SCode == "ref")
for (ii in levels(ara.df$Batch))
  lines(ara.df$SeqNr[ara.df$Batch == ii],
        predict(BLines2, newdata = ara.df[ara.df$Batch == ii, ]),
        col = 2)

correctfun <- function(x, seqnr, batch) {
  huhn.df <- data.frame(x = x, seqnr = seqnr, batch = batch)
  blines <- lm(x ~ seqnr * batch, data = huhn.df)
  huhn.df$qcpred <- predict(blines, newdata = huhn.df)
  huhn.df$x - huhn.df$qcpred
}

nna.threshold <- 50
x.idx <- apply(arabidopsis, 2, function(x) sum(is.na(x)) < 50)
correctedX <- apply(arabidopsis[, x.idx], 2, correctfun,
                    seqnr = arabidopsis.Y$SeqNr,
                    batch = arabidopsis.Y$Batch)

XsampleC <- correctedX
for (i in 1:ncol(XsampleC))
  XsampleC[is.na(XsampleC[, i]), i] <- mean(XsampleC[, i], na.rm = TRUE)
Xc.PCA <- PCA(scale(XsampleC))
bhatta.corr <- aveBhatta(scores(Xc.PCA, 2), arabidopsis.Y$Batch)

X.PCA2 <- PCA(scale(Xsample[, x.idx]))
bhatta.origsub <- aveBhatta(scores(X.PCA2, 2), arabidopsis.Y$Batch)

scoreplot(Xc.PCA, col = cols, pch = pchs)
legend("topleft", legend = levels(arabidopsis.Y$Batch),
       ncol = 2, pch = 1:nlevels(arabidopsis.Y$Batch), col = lcols)
title(main = paste("Average between-batch distance:",
                   round(bhatta.corr, 3)))

X <- arabidopsis[, x.idx]
na.idx <- is.na(X)
X[na.idx] <- min(X, na.rm = TRUE)

require(RUVSeq, quiet = TRUE)
idx <- which(arabidopsis.Y$SCode == "ref")
replicates.ind <-
  matrix(-1, nrow(X) - length(idx) + 1, length(idx))
replicates.ind[1, ] <- idx
replicates.ind[-1, 1] <- (1:nrow(X))[-idx]
X.RUVcorrected <-
  t(RUVs(x = t(X), cIdx = 1:ncol(X), k = 3,
         scIdx = replicates.ind, isLog = TRUE)$normalizedCounts)
X.RUVcorrected[na.idx] <- NA

XR <- X.RUVcorrected
for (i in 1:ncol(XR))
  XR[is.na(XR[, i]), i] <- mean(XR[, i], na.rm = TRUE)
XRUV.PCA <- PCA(scale(XR))
bhatta.RUV <- aveBhatta(scores(XRUV.PCA, 2), arabidopsis.Y$Batch)

data(bdata)
par(mar = c(0, 3, 0, 0), mfrow = c(1, 2))
persp(bdata$d1, phi = 20, theta = 34, expand = .5,
      xlab = "Time", ylab = "Wavelength")
par(mar = c(7.1, 4.5, 4.1, 4.5))
matplot(cbind(c(bdata$sp1), c(bdata$sp2)), type = "l", lty = 1,
        xlab = "Wavelength number", ylab = "Intensity")

efa <- function(x, ncomp) {
  nx <- nrow(x)
  Tos <- Fros <- matrix(0, nx, ncomp)
  for (i in 3:nx)
    Tos[i, ] <- svd(scale(x[1:i, ], scale = FALSE))$d[1:ncomp]
  for (i in (nx-2):1)
    Fros[i, ] <- svd(scale(x[i:nx, ], scale = FALSE))$d[1:ncomp]

  Combos <- array(c(Tos, Fros[, ncomp:1]), c(nx, ncomp, 2))
  list(forward = Tos, backward = Fros,
       pure.comp = apply(Combos, c(1, 2), min))
}

par(opar)
par(mfrow = c(1, 2))
X <- bdata$d1
X.efa <- efa(X, 3)
matplot(X.efa$forward, type = "l", ylab = "Singular values", lty=1)
matplot(X.efa$backward, type = "l", ylab = "Singular values", lty=1)

par(opar)
matplot(X.efa$pure.comp, type = "l", ylab = "", lty = 1)
legend("topright", legend = paste("Comp", 1:3),
       lty = 1, col = 1:3, bty="n")

require(alsace, quiet = TRUE)
X.opa <- opa(X, 3)

matplot(X.opa, type = "l", col = 1:3, lty = 1,
        ylab = "response", xlab = "wavelength number")
legend("topleft", legend = paste("Comp", 1:3), col = 1:3,
       lty = 1, bty = "n")

normS <- function(S) {
  sweep(S,
        1,
        apply(S, 1, function(x) sqrt(sum(x^2))),
        FUN = "/")
}

getS <- function(data, C) {
  normS(ginv(C) %*% data)
}
getC <- function(data, S) {
  data %*% ginv(S)
}

mcr <- function(x, init, what = c("row", "col"),
                convergence = 1e-8, maxit = 50) {
  what <- match.arg(what)
  if (what == "col") {
    CX <- init
    SX <- getS(x, CX)
  } else {
    SX <- normS(init)
    CX <- getC(x, SX)
  }

  rms <- rep(NA, maxit + 1)
  rms[1] <- sqrt(mean((x - CX %*% SX)^2))

  for (i in 1:maxit) {
    CX <- getC(x, SX)
    SX <- getS(x, CX)

    resids <- x - CX %*% SX
    rms[i+1] <- sqrt(mean(resids^2))
    if ((rms[i] - rms[i+1]) < convergence) break;
  }

  list(C = CX, S = SX, resids = resids, rms = rms[!is.na(rms)])
}

par(mfrow = c(1, 2))
X.mcr.efa <- mcr(X, X.efa$pure.comp, what = "col")
matplot(X.mcr.efa$C, type = "n", 
        main = "Concentration profiles", ylab = "Concentration")
abline(h = 0, col = "lightgray")
matlines(X.efa$pure.comp * 5, type = "l", lty = 2, col = 1:3)
matlines(X.mcr.efa$C, type = "l", col = 1:3, lty = 1, lwd = 2)
legend("topright", legend = paste("Comp", 1:3),
       col = 1:3, lty = 1, bty = "n")
matplot(t(X.mcr.efa$S), col = 1:3, type = "l", lty = 1,
        main = "Pure spectra", ylab = "Intensity")
abline(h = 0, col = "lightgray")
legend("topright", legend = paste("Comp", 1:3), lty = 1,
       bty = "n", col = 1:3) 

X.mcr.opa <- mcr(X, t(X.opa), what = "row")

X.mcr.efa$rms
X.mcr.opa$rms

X.als.efa <-  als(CList = list(X.efa$pure.comp),
                  PsiList = list(X), S = matrix(0, 73, 3), 
                  nonnegS = TRUE, nonnegC = TRUE,
                  normS = 0.5, uniC = TRUE)

X.als.opa <- als(CList = list(matrix(0, 40, 3)),
                 PsiList = list(X), S = X.opa,
                 nonnegS = TRUE, nonnegC = TRUE,
                 optS1st = FALSE, normS = 0.5, uniC = TRUE)

par(mfrow = c(2, 2))
matplot(X.als.efa$S, type = "n", main = "Pure spectra (EFA)",
        ylab = "Intensity")
abline(h = 0, col = "lightgray")
matlines(X.als.efa$S, type = "l", lty = 1, col = 1:3)
legend("topright", legend = paste("Comp", 1:3), lty = 1,
       bty = "n", col = 1:3) 
matplot(X.als.efa$CList[[1]], type = "n", 
        main = "Concentration profiles (EFA)",
        ylab = "Concentration")
abline(h = 0, col = "lightgray")
matlines(X.als.efa$CList[[1]], lty = 1, col = 1:3)

matplot(X.als.opa$S, type = "n", main = "Pure spectra (OPA)",
        ylab = "Intensity")
abline(h = 0, col = "lightgray")
matlines(X.als.opa$S, type = "l", lty = 1, col = 1:3)
legend("topright", legend = paste("Comp", 1:3), lty = 1,
       bty = "n", col = 1:3) 
matplot(X.als.opa$CList[[1]], type = "n", 
        main = "Concentration profiles (OPA)",
        ylab = "Concentration")
abline(h = 0, col = "lightgray")
matlines(X.als.opa$CList[[1]], lty = 1, col = 1:3)

C0 <- matrix(0, 40, 3)
X2.als.opa <-
  als(CList = list(C0, C0),
      PsiList = list(bdata$d1, bdata$d2),
      S = X.opa, normS = 0.5,
      nonnegS = TRUE, nonnegC = TRUE,
      optS1st = FALSE, uniC = TRUE)

cor(X.als.opa$S, cbind(c(bdata$sp1), c(bdata$sp2)))
cor(X2.als.opa$S, cbind(c(bdata$sp1), c(bdata$sp2)))

