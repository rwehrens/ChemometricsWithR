## Demo featuring all R code from Chapter 9 in the book "Chemometrics
## with R" by Ron Wehrens. The only additions to the book chapter code
## are loading the required libraries and data, and occasionally
## closing a graphics window.

graphics.off()
opar <- par(no.readonly = TRUE)

require(pls)
data(gasoline, package = "pls")
wavelengths <- seq(900, 1700, by = 2)
gas.odd <- seq(1, nrow(gasoline$NIR), by = 2)
gas.even <- seq(2, nrow(gasoline$NIR), by = 2)

gasoline.mscpcr <- pcr(octane ~ msc(NIR), data = gasoline, 
                       ncomp = 4)
gasoline.mscpcr.cv <- crossval(gasoline.mscpcr, length.seg = 1)
RMSEP(gasoline.mscpcr.cv, estimate = "CV")

require(class)
data(wines, package = "kohonen")
wine.classes <- as.integer(vintages)
wines.sc <- scale(wines)
wines.odd <- seq(1, nrow(wines), by = 2)
wines.even <- seq(2, nrow(wines), by = 2)
twowines <- wines[vintages != "Barolo",]
twovintages <- factor(vintages[vintages != "Barolo"])

X <- wines[vintages != "Barolo", ]
vint <- factor(vintages[vintages != "Barolo"])
kvalues <- 1:12
set.seed(17)
ktabs <- lapply(kvalues,
                function(i) {
                  kpred <- knn.cv(X, vint, k = i)
                  table(vint, kpred)
                })

TPrates <- sapply(ktabs, function(x) x[1, 1]/sum(x[, 1]))
FPrates <- sapply(ktabs, function(x) 1 - x[2, 2]/sum(x[2, ]))
plot(FPrates, TPrates, type = "b",
     xlim = c(.15, .45), ylim = c(.5, .75),
     xlab = "FP rate", ylab = "TP rate")
text(FPrates, TPrates, 1:12, pos = 4)

gasoline.pcr <- pcr(octane ~ ., data = gasoline,
                    validation = "LOO", ncomp = 1)
RMSEP(gasoline.pcr, estimate = "CV")

gasoline.pcr2 <- pcr(octane ~ ., data = gasoline, ncomp = 1)
X <- gasoline.pcr2$scores
HatM <- X %*% solve(crossprod(X), t(X))
sqrt(mean((gasoline.pcr2$residuals/(1 - diag(HatM)))^2))

sqrt(mean((gasoline.pcr2$residuals/(1 - mean(diag(HatM))))^2))

gasoline.pcr <- pcr(octane ~ ., data = gasoline,
                    validation = "CV", ncomp = 4,
                    segment.type = "interleaved")
RMSEP(gasoline.pcr, estimate = "CV")

gasoline.pls <- plsr(octane ~ ., data = gasoline,
                     validation = "LOO", ncomp = 2,
                     jackknife = TRUE)
n <- length(gasoline$octane)
b.oob <- gasoline.pls$validation$coefficients[, , 2, ]
bias.est <- (n-1) * (rowMeans(b.oob) - coef(gasoline.pls))
plot(wavelengths, bias.est, xlab = "wavelength", ylab = "bias",
     type = "h", main = "Jackknife bias estimates",
     col = "gray")
var.est <- var.jack(gasoline.pls)
lines(wavelengths, var.est, col = "red")

set.seed(7)
B <- 500
ngas <- nrow(gasoline)
boot.indices <- 
  matrix(sample(1:ngas, ngas * B, replace = TRUE), ncol = B)
sort(boot.indices[, 1])[1:20]

set.seed(7)
npc <- 5
predictions <- array(NA, c(ngas, npc, B))
for (i in 1:B) {
  gas.bootpcr <- pcr(octane ~ ., data = gasoline,
                     ncomp = npc, subset = boot.indices[, i])
  oobs <- (1:ngas)[-boot.indices[, i]]
  predictions[oobs, , i] <- 
    predict(gas.bootpcr,
            newdata = gasoline$NIR[oobs, ])[, 1, ]
}

diffs <- sweep(predictions, 1, gasoline$octane)
sqerrors <- apply(diffs^2, c(1, 2), mean, na.rm = TRUE)
sqrt(colMeans(sqerrors))

gas.pcr <- pcr(octane ~ ., data = gasoline, ncomp = npc)
RMSEP(gas.pcr, intercept = FALSE)
error.632 <- .368 * colMeans(gas.pcr$residuals^2) +
  .632 * colMeans(sqerrors)
sqrt(error.632)

gas.pcr.cv <- pcr(octane ~ ., data = gasoline, ncomp = npc,
                  validation = "CV")
gas.pcr.loo <- pcr(octane ~ ., data = gasoline, ncomp = npc,
                   validation = "LOO")
bp <- barplot(sqrt(error.632),
              ylim = c(0, 1.6), col = "peachpuff")
lines(bp, sqrt(c(gas.pcr.cv$validation$PRESS) / ngas), 
      col = 2, lwd = 2)
lines(bp, sqrt(c(gas.pcr.loo$validation$PRESS) / ngas), 
      col = 4, lty = 2, lwd = 2)
legend("topright", lty = 1:2, col = c(2, 4), lwd = 2,
       legend = c("CV", "LOO"))

require(boot)
set.seed(7)
gas.pcr.boot632 <- 
  boot(gasoline, 
       function(x, ind) {
         mod <- pcr(octane ~ ., data = x,
                    subset = ind, ncomp = 4)
         gasoline$octane - 
           predict(mod, newdata = gasoline$NIR, ncomp = 4)},
       R = 499)

dim(boot.array(gas.pcr.boot632))
boot.array(gas.pcr.boot632)[1, 1:10]

in.bag <- boot.array(gas.pcr.boot632)
oob.error <- mean((gas.pcr.boot632$t^2)[in.bag == 0])
app.error <- MSEP(pcr(octane ~ ., data = gasoline, ncomp = 4),
                  ncomp = 4, intercept = FALSE)
sqrt(.368 * c(app.error$val) + .632 * oob.error)

B <- 1000
ngas <- nrow(gasoline)
boot.indices <- 
  matrix(sample(1:ngas, ngas * B, replace = TRUE), ncol = B)
npc <- 4
gas.pcr <- pcr(octane ~ ., data = gasoline, ncomp = npc)
coefs <- matrix(0, ncol(gasoline$NIR), B)
for (i in 1:B) {
  gas.bootpcr <- pcr(octane ~ ., data = gasoline,
                     ncomp = npc, subset = boot.indices[, i])
  coefs[, i] <- c(coef(gas.bootpcr))
}

matplot(wavelengths, coefs, type = "n",
        ylab = "Coefficients", xlab = "Wavelength (nm)")
abline(h = 0, col = "gray")
polygon(c(wavelengths, rev(wavelengths)),
        c(apply(coefs, 1, max), rev(apply(coefs, 1, min))),
        col = "steelblue", border = NA)

coef.stats <- cbind(apply(coefs, 1, quantile, .025),
                    apply(coefs, 1, quantile, .975))
matplot(wavelengths, coef.stats, type = "n",
        xlab = "Wavelength (nm)",
        ylab = "Regression coefficient")
abline(h = 0, col = "gray")
polygon(c(wavelengths, rev(wavelengths)),
        c(coef.stats[, 1], rev(coef.stats[, 2])),
        col = "pink", border = NA)
lines(wavelengths, c(coef(gas.pcr)))

set.seed(17)
gas.pcr.bootCI <-
  boot(gasoline,
       function(x, ind) {
         c(coef(pcr(octane ~ ., data = x,
                    ncomp = npc, subset = ind)))},
       R = 999)
dim(gas.pcr.bootCI$t)

smallest <- which.min(gas.pcr.bootCI$t0)
plot(gas.pcr.bootCI, index = smallest)

boot.ci(gas.pcr.bootCI, index = smallest, type = c("perc", "bca"))

all.cis <- lapply(1:length(wavelengths),
                  function(ii)
                    boot.ci(gas.pcr.bootCI, index = ii,
                            type = c("perc", "bca")))
perc.cis <- t(sapply(all.cis, "[[", "percent"))
bca.cis <- t(sapply(all.cis, "[[", "bca"))
nperc <- sum(perc.cis[, 4] > 0) + sum(perc.cis[, 5] < 0)
nbca <- sum(bca.cis[, 4] > 0) + sum(bca.cis[, 5] < 0)

set.seed(7)
gas.BIG.bootCI <-
  boot(gasoline,
       function(x, ind) {
         c(coef(pcr(octane ~ ., data = x,
                    ncomp = 20, subset = ind)))},
       R = 999)
big.cis <- lapply(1:length(wavelengths),
                  function(ii)
                    boot.ci(gas.BIG.bootCI, index = ii,
                            type = c("perc", "bca")))
PERC.cis <- t(sapply(big.cis, "[[", "percent"))
BCA.cis <- t(sapply(big.cis, "[[", "bca"))
NPERC <- sum(PERC.cis[, 4] > 0) + sum(PERC.cis[, 5] < 0)
NBCA <- sum(BCA.cis[, 4] > 0) + sum(BCA.cis[, 5] < 0)

require(ipred)
set.seed(77)
(gasoline.bagging <- ipredbagg(gasoline$octane[gas.odd],
                               gasoline$NIR[gas.odd, ],
                               coob = TRUE))

gs.baggpreds <-
  predict(gasoline.bagging, gasoline$NIR[gas.even, ])
resids <- gs.baggpreds - gasoline$octane[gas.even]
sqrt(mean(resids^2))

data(Prostate2000Raw, package = "ChemometricsWithR")                            
prostate <- rowsum(t(Prostate2000Raw$intensity),                                
                     group = rep(1:327, each = 2),                              
                     reorder = FALSE) / 2                                       
prostate.type <- Prostate2000Raw$type[seq(1, 654, by = 2)]                      

require(kohonen, quiet = TRUE) # for classvec2classmat
prost <- prostate[prostate.type != "bph", 1:1000]
prost.type <- factor(prostate.type[prostate.type != "bph"])
prost.df <- data.frame(type = prost.type, prost = prost)
prost.odd <- seq(1, length(prost.type), by = 2)
prost.even <- seq(2, length(prost.type), by = 2)
prostate.clmat <- classvec2classmat(prostate.type)
prostate.df <- data.frame(class = I(prostate.clmat),
                          msdata = I(prostate[, 1:1000]))
prostate.odd <- seq(1, nrow(prostate.df), by = 2)
prostate.even <- seq(2, nrow(prostate.df), by = 2)

set.seed(78)
prost.bagging <- bagging(type ~ ., data = prost.df, 
                         subset = prost.odd)
prost.baggingpred <- predict(prost.bagging, 
                             newdata = prost.df[prost.even, ])
table(prost.type[prost.even], prost.baggingpred)

require(randomForest, quiet = TRUE)
wines.df <- data.frame(vint = vintages, wines)
(wines.rf <- randomForest(vint ~ ., subset = wines.odd,
                          data = wines.df))

wines.rf.predict <-
  predict(wines.rf, newdata = wines.df[wines.even, ])
table(wines.rf.predict, vintages[wines.even])

wines.rf <- randomForest(vint ~ ., data = wines.df,
                         importance = TRUE)
varImpPlot(wines.rf)

gasoline.rf <- randomForest(gasoline$NIR[gas.odd, ], 
                            gasoline$octane[gas.odd],
                            importance = TRUE,
                            xtest = gasoline$NIR[gas.even, ],
                            ytest = gasoline$octane[gas.even])

par(mfrow = c(1, 2))
pl.range <- c(83, 90)
plot(gasoline$octane[gas.odd], gasoline.rf$predicted,
     main = "Training: OOB prediction", xlab = "True",
     ylab = "Predicted", xlim = pl.range, ylim = pl.range)
abline(0, 1, col = "gray")
plot(gasoline$octane[gas.even], gasoline.rf$test$predicted,
     main = "Test set prediction", xlab = "True",
     ylab = "Predicted", xlim = pl.range, ylim = pl.range)
abline(0, 1, col = "gray")

resids <- gasoline.rf$test$predicted - gasoline$octane[gas.even]
sqrt(mean(resids^2))

par(opar)
rf.imps <- importance(gasoline.rf)
plot(wavelengths, rf.imps[, 1] / max(rf.imps[, 1]),
     type = "l", xlab = "Wavelength (nm)",
     ylab = "Importance", col = "gray")
lines(wavelengths, rf.imps[, 2] / max(rf.imps[, 2]), col = 2)
legend("topright", legend = c("Error decrease", "Gini index"),
       col = c("gray", "red"), lty = 1)

prost.rf <-
  randomForest(x = prost[prost.odd, ],
               y = prost.type[prost.odd],
               x.test = prost[prost.even, ],
               y.test = prost.type[prost.even])
prost.rfpred <- predict(prost.rf, newdata = prost[prost.even, ])
table(prost.type[prost.even], prost.rfpred)

require(ada, quiet = TRUE)
prost.ada <- ada(type ~ ., data = prost.df, subset = prost.odd)
prost.adapred <-
  predict(prost.ada, newdata = prost.df[prost.even, ])
table(prost.type[prost.even], prost.adapred)

prost.ada <- addtest(prost.ada,
                     prost.df[prost.even, ],
                     prost.type[prost.even])
plot(prost.ada, test = TRUE)

