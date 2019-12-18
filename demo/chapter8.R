## Demo featuring all R code from Chapter 8 in the book "Chemometrics
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

X <- gasoline$NIR[, 100*(1:4)]
Y <- gasoline$octane
Xtr <- cbind(1, X[gas.odd, ])
Ytr <- Y[gas.odd]
Bs <- t(solve(crossprod(Xtr), t(Xtr)) %*% Ytr)
Bs

gas.df <- data.frame(octane = gasoline$octane,
                     V100 = gasoline$NIR[, 100],
                     V200 = gasoline$NIR[, 200],
                     V300 = gasoline$NIR[, 300],
                     V400 = gasoline$NIR[, 400])
Blm <- lm(octane ~ ., data = gas.df)
summary(Blm)

X <- scale(gasoline$NIR, scale = FALSE,
           center = colMeans(gasoline$NIR[gas.odd, ]))

Xgas.odd.svd <- svd(X[gas.odd, ])
Xgas.odd.scores <- Xgas.odd.svd$u %*% diag(Xgas.odd.svd$d)
gas.odd.pcr <- 
  lm(gasoline$octane[gas.odd] ~ I(Xgas.odd.scores[, 1:5]) - 1)

gas.odd.coefs <- coef(gas.odd.pcr) %*% t(Xgas.odd.svd$v[, 1:5])

gasoline.pcr <- pcr(octane ~ ., data = gasoline, 
                    subset = gas.odd, ncomp = 5)
all.equal(c(coef(gasoline.pcr)), c(gas.odd.coefs))
plot(wavelengths, coef(gasoline.pcr), type = "l",
     xlab = "Wavelength (nm)", ylab = "Regression coefficient")
abline(h = 0, col = "gray")

summary(gasoline.pcr)

RMSEP(gasoline.pcr, estimate = "train", intercept = FALSE)

RMSEP(gasoline.pcr, estimate = "train", comp = 4)

par(mfrow = c(1, 2))
gasoline.pcr <- pcr(octane ~ ., data = gasoline, subset = gas.odd,
                    validation = "LOO", ncomp = 10)
plot(gasoline.pcr, "validation", estimate = "CV")
par(pty = "s")
plot(gasoline.pcr, "prediction", ncomp = 4)
abline(0, 1, col = "gray")

RMSEP(gasoline.pcr, estimate = "CV")

sqrt(gasoline.pcr$validation$PRESS / nrow(Xtr))

gasoline.pcr.pred <- predict(gasoline.pcr, ncomp = 4,
                             newdata = gasoline[gas.even, ])
rms(gasoline$octane[gas.even], gasoline.pcr.pred)

RMSEP(gasoline.pcr, ncomp = 4, newdata = gasoline[gas.even, ],
      intercept = FALSE)

par(mfrow = c(1, 2))
nc.1s <-
  selectNcomp(gasoline.pcr, method = "onesigma", plot = TRUE)
nc.rand <-
  selectNcomp(gasoline.pcr, method = "randomization", plot = TRUE)

gasoline.pls <- plsr(octane ~ ., data = gasoline, 
                     subset = gas.odd, ncomp = 5)
summary(gasoline.pls)

gasoline.pls <- plsr(octane ~ ., data = gasoline, subset = gas.odd,
                     validation = "LOO", ncomp = 10)
par(mfrow = c(1, 2))
plot(gasoline.pls, "validation", estimate = "CV")
par(pty = "s")
plot(gasoline.pls, "prediction", ncomp = 3)
abline(0, 1, col = "gray")

RMSEP(gasoline.pls, ncomp = 3, newdata = gasoline[gas.even, ],
      intercept = FALSE)

cor(gasoline.pls$loadings[, 1:3])

cor(gasoline.pls$scores[, 1:3])

par(opar)
plot(gasoline.pls, "loading", comps = 1:3, legendpos = "topleft",
     lty = 1, col = 1:3)

par(mfrow = c(1, 2))
plot(scores(gasoline.pls)[, 1], Yscores(gasoline.pls)[, 1],
     xlab = "X scores", ylab = "Y scores", main = "LV 1")
abline(h = 0, v = 0, col = "gray")
plot(scores(gasoline.pls)[, 2], Yscores(gasoline.pls)[, 2],
     xlab = "X scores", ylab = "Y scores", main = "LV 2")
abline(h = 0, v = 0, col = "gray")

require(MASS)
gasoline.ridge <- 
  lm.ridge(octane ~ NIR, data = gasoline, subset = gas.odd,
           lambda = seq(0.001, 0.1, by = 0.01))
select(gasoline.ridge)

par(mfrow = c(1, 2))

plot(-5:5, -5:5, ylim = c(-.2, 4), xlim = c(-3, 3), axes = FALSE, type = "n",
     ylab = "", xlab = "")
mtext("V(e)", 2, line = 1)
mtext("e", 1, line = 2)
mtext(c(expression(-epsilon), expression(epsilon)), at = c(-1, 1),
      line = .5, side = 1)
box()
abline(h = 0, v = 0, col = "gray")
abline(v = c(-1, 1), col = "gray", lty = 2)
segments(c(-6, -1, 1),
         c(6, 0, 0),
         c(-1, 1, 6),
         c(0, 0, 6))

plot(-5:5, -5:5, ylim = c(-.2, 4), xlim = c(-3, 3), axes = FALSE, type = "n",
     ylab = "", xlab = "")
mtext("V(e)", 2, line = 1)
mtext("e", 1, line = 2)
mtext(c(expression(-epsilon), expression(epsilon)), at = c(-1, 1),
      line = .5, side = 1)
box()
abline(h = 0, v = 0, col = "gray")
abline(v = c(-1, 1), col = "gray", lty = 2)

segments(-6, 6, -1, .5)
segments(6, 6, 1, .5)
x <- seq(-1, 1, length = 100)
y <- .5 * x^2
lines(x, y)

require(e1071)
set.seed(7)
gasoline.svm <- svm(octane ~ ., data = gasoline, 
                    subset = gas.odd, cross = 10)
(gasoline.svmsum <- summary(gasoline.svm))

par(mfrow = c(1, 2))
plot(gasoline$octane[gas.odd], predict(gasoline.svm),
     main = "Training set", xlab = "Octane number (true)", 
     ylab = "Octane number (predicted)")
abline(0, 1)
plot(gasoline$octane[gas.even], 
     predict(gasoline.svm, new = gasoline[gas.even, ]),
     main = "Test set", xlab = "Octane number (true)", 
     ylab = "Octane number (predicted)")
abline(0, 1)

rms(gasoline$octane[gas.even], 
    predict(gasoline.svm, new = gasoline[gas.even, ]))

gasoline.svmlin <- svm(octane ~ ., data = gasoline, 
                       subset = gas.odd, kernel = "linear")
rms(gasoline$octane[gas.even], 
    predict(gasoline.svmlin, new = gasoline[gas.even, ]))

require(nnet)
X <- scale(gasoline$NIR, scale = FALSE,
           center = colMeans(gasoline$NIR[gas.odd, ]))
X.odd.svd <- svd(X[gas.odd, ])
X.odd.scores <- X.odd.svd$u %*% diag(X.odd.svd$d)
X.even.scores <- X[gas.even, ] %*% X.odd.svd$v
gas.nnet <- nnet(X.odd.scores[, 1:5],
                 matrix(gasoline$octane[gas.odd], ncol = 1),
                 size = 5, linout = TRUE)

gas.nnet.pred <- predict(gas.nnet, X.even.scores)
rms(gas.nnet.pred, gasoline$octane[gas.even])

par(opar)
set.seed(7)
rms.values <- rep(0, 100)
for (i in 1:100) {
  gas.nnet <- nnet(Xgas.odd.scores[, 1:5],
                   matrix(gasoline$octane[gas.odd], ncol = 1),
                   size = 5, linout = TRUE)
  gas.nnet.pred <- predict(gas.nnet, X.even.scores)
  gas.nnet.res <- gas.nnet.pred - gasoline$octane[gas.even]
  rms.values[i] <- sqrt(mean(gas.nnet.res^2))
}
hist(rms.values, breaks = 50, xlim = c(0, 1.8),
     main = "Repeated ANN training", , xlab = "RMS")

## Wine data
require(kohonen) # for classvec2classmat, and the data
data(wines, package = "kohonen")
wine.classes <- as.integer(vintages)
wines.sc <- scale(wines)
wines.odd <- seq(1, nrow(wines), by = 2)
wines.even <- seq(2, nrow(wines), by = 2)
twowines <- wines[vintages != "Barolo",]
twovintages <- factor(vintages[vintages != "Barolo"])

C <- classvec2classmat(vintages[c(1:25, 61:85)])
X <- wines[c(1:25, 61:85), c(7, 13)]

C <- classvec2classmat(vintages[c(1:25, 61:85)])
X <- wines[c(1:25, 61:85), c(7, 13)]
wines.lm <- lm(C ~ X)
wines.lm.predict <- classmat2classvec(predict(wines.lm))
table(vintages[c(1:25, 61:85)], wines.lm.predict)

wines.lda <- lda(factor(vintages[c(1:25, 61:85)]) ~ X)
table(vintages[c(1:25, 61:85)], predict(wines.lda)$class)

## Prostate data
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

prost.tcp <- tcrossprod(scale(prost))
prost.svd <- svd(prost.tcp)
prost.scores <- prost.svd$u %*% diag(sqrt(prost.svd$d))

pairs(prost.scores[, 1:5], labels = paste("PC", 1:5),
      pch = as.integer(prost.type), col = as.integer(prost.type))

prost.pcda5 <- lda(prost.type ~ prost.scores[, 1:5], CV = TRUE) # INCORRECT
(prost.ldaresult <- table(prost.type, prost.pcda5$class))

par(mfrow = c(1, 2))
set.seed(7)
prost.df2 <- data.frame(class = as.integer(prost.type),
                       msdata = I(prost))
prost.pcr <- pcr(class ~ msdata, ncomp = 16,
                 data = prost.df2, subset = prost.odd,
                 validation = "CV", scale = TRUE)
validationplot(prost.pcr, estimate = "CV")
prost.cv.cl <- round(prost.pcr$validation$pred[, 1, ])
prost.cv.err <- apply(prost.cv.cl, 2, 
                      err.rate, prost.df2$class[prost.odd])
prost.tst <- predict(prost.pcr, newdata = prost.df2[prost.even, ])
prost.tst.cl <- round(prost.tst[, 1, ])
prost.tst.err <- apply(prost.tst.cl, 2, 
                       err.rate, prost.df2$class[prost.even])

matplot(cbind(prost.cv.err, prost.tst.err),
        lty = 1:2, col = 2:1, type = "l",
        xlab = "# PCs", ylab = "Misclassif. rate")
legend("topright", lty = 1:2, col = 2:1, bty = "n",
       legend = c("Training set CV", "Test set prediction"))

set.seed(17)
prostate.pcr <- pcr(class ~ msdata, ncomp = 16,
                    data = prostate.df,
                    subset = prostate.odd,
                    validation = "CV", scale = TRUE)

pcr.predictions.loo <-
  sapply(1:16,
         function(i, arr) classmat2classvec(arr[, , i]),
         prostate.pcr$validation$pred)
pcr.loo.err <- apply(pcr.predictions.loo, 2, err.rate, 
                     prostate.type[prostate.odd])

prostate.pcrpred <- 
  predict(prostate.pcr, new = prostate.df[prost.even, ])
predictions.pcrtest <-
  sapply(1:16,
         function(i, arr) classmat2classvec(arr[, , i]),
         prostate.pcrpred)

table(prostate.type[prost.even], predictions.pcrtest[, 7])

set.seed(17)
prostate.pls <- plsr(class ~ msdata, data = prostate.df,
                     subset = prostate.odd, scale = TRUE,
                     ncomp = 16, validation = "CV") 

pls.predictions.loo <-
  sapply(1:16,
         function(i, arr) classmat2classvec(arr[, , i]),
         prostate.pls$validation$pred)
pls.loo.err <- apply(pls.predictions.loo, 2, err.rate, 
                     prostate.type[prostate.odd])
prostate.plspred <- 
  predict(prostate.pls, new = prostate.df[prostate.even, ])
predictions.plstest <-
  sapply(1:16,
         function(i, arr) classmat2classvec(arr[, , i]),
         prostate.plspred)

table(prostate.type[prostate.even], predictions.plstest[, 4])

par(mfrow = c(1, 2))
matplot(cbind(pcr.loo.err,
              apply(predictions.pcrtest, 2, err.rate, 
                    prostate.type[prost.even])),
        type = "l", lty = 1:2, col = 2:1,
        main = "PCDA",
        xlab = "# PCs", ylab = "Misclassif. rate")
legend("topright", lty = 1:2, col = 2:1, bty = "n",
       legend = c("Training set CV", "Test set prediction"))
pls.predictions.loo <-
  sapply(1:16,
         function(i, arr) classmat2classvec(arr[, , i]),
         prostate.pls$validation$pred)
pls.loo.err <- apply(pls.predictions.loo, 2, err.rate, 
                     prostate.type[prostate.odd])
prostate.plspred <- 
  predict(prostate.pls, new = prostate.df[prostate.even, ])
predictions.plstest <-
  sapply(1:16,
         function(i, arr) classmat2classvec(arr[, , i]),
         prostate.plspred)
matplot(cbind(pls.loo.err,
              apply(predictions.plstest, 2, err.rate, 
                    prostate.type[prostate.even])),
        type = "l", lty = 1:2, col = 2:1,
        main = "PLSDA", xlab = "# PCs", ylab = "Misclassif. rate")
legend("topright", lty = 1:2, col = 2:1, bty = "n",
       legend = c("Training set CV", "Test set prediction"))

prostate.ldapls <- lda(scores(prostate.pls)[, 1:6], 
                       prostate.type[prostate.odd])

tst.scores <- predict(prostate.pls, ncomp = 1:6,
                      newdata = prostate.df[prostate.even, ],
                      type = "scores")
table(prostate.type[prostate.even],
      predict(prostate.ldapls, new = tst.scores)$class)

nvar <- 2000
nobj <- 40
RandX <- matrix(rnorm(nobj*nvar), nrow = nobj)
RandY <- sample(c(0, 1), nobj, replace = TRUE)

Rand.pcr <- pcr(RandY ~ RandX, ncomp = 2)
Rand.ldapcr <- lda(RandY ~ scores(Rand.pcr), CV = TRUE)
table(RandY, Rand.ldapcr$class)

Rand.pls <- plsr(RandY ~ RandX, ncomp = 2)
Rand.ldapls <- lda(RandY ~ scores(Rand.pls), CV = TRUE)
table(RandY, Rand.ldapls$class)

