## Demo featuring all R code from Chapter 4 in the book "Chemometrics
## with R" by Ron Wehrens. The only additions to the book chapter code
## are loading the required libraries and data, and occasionally
## closing a graphics window.
graphics.off()
opar <- par(no.readonly = TRUE)

data(wines, package = "kohonen")
wine.classes <- as.integer(vintages)
wines.sc <- scale(wines)
wines.svd <- svd(wines.sc)
wines.scores <- wines.svd$u %*% diag(wines.svd$d)
wines.loadings <- wines.svd$v

wines.vars <- wines.svd$d^2 / (nrow(wines) - 1)
wines.totalvar <- sum(wines.vars)
wines.relvars <- wines.vars / wines.totalvar
variances <- 100 * round(wines.relvars, digits = 3)
variances[1:5]

par(mfrow = c(1, 2))
plot(wines.scores[, 1:2], type = "n",
     xlab = paste("PC 1 (", variances[1], "%)", sep = ""),
     ylab = paste("PC 2 (", variances[2], "%)", sep = ""))
abline(h = 0, v = 0, col = "gray")
points(wines.scores[, 1:2], pch = wine.classes, col = wine.classes)
plot(wines.loadings[, 1] * 1.2, wines.loadings[, 2], type = "n",
     xlab = paste("PC 1 (", variances[1], "%)", sep = ""),
     ylab = paste("PC 2 (", variances[2], "%)", sep = ""))
abline(h = 0, v = 0, col = "gray", lty = 2)
arrows(0, 0, wines.loadings[, 1], wines.loadings[, 2],
     col = "blue", length = .15, angle = 20)
text(wines.loadings[, 1:2], labels = colnames(wines))

nr <- 100
nc <- 100000
X <- matrix(rnorm(nr*nc), nrow = nr)
system.time({
  X.svd <- svd(X)
  X.scores <- X.svd$u %*% diag(X.svd$d)
  X.variances <- X.svd$d^2 / (nrow(X) - 1)
  X.loadings <- X.svd$v
})

system.time({
  X2.tcp <- tcrossprod(X)
  X2.svd <- svd(X2.tcp)
  X2.scores <- X2.svd$u %*% diag(sqrt(X2.svd$d))
  X2.variances <- X2.svd$d / (nrow(X) - 1)
  X2.loadings <- solve(X2.scores, X)
})

par(mfrow = c(2, 2))
barplot(wines.vars[1:10], main = "Variances",
        names.arg = paste("PC", 1:10))
barplot(log(wines.vars[1:10]), main = "log(Variances)",
        names.arg = paste("PC", 1:10))
barplot(wines.relvars[1:10], main = "Relative variances",
        names.arg = paste("PC", 1:10))
barplot(cumsum(100 * wines.relvars[1:10]), 
        main = "Cumulative variances (%)",
        names.arg = paste("PC", 1:10), ylim = c(0, 100))

llambdas <- log(wines.vars)
CIwidth <- qnorm(.975) * sqrt(2 / (nrow(wines) - 1))
CIs <- cbind(exp(llambdas - CIwidth),
             wines.vars,
             exp(llambdas + CIwidth))
colnames(CIs) <- c("CI 0.025", "  Estimate", "  CI 0.975")
CIs[1:4, ]

small.ones <- wines.vars[11:13]
n <- nrow(wines)
nsmall <- length(small.ones)
geo.mean <- prod(small.ones)^{1/nsmall}
mychisq <- (n - 1) * nsmall * log(mean(small.ones) / geo.mean)
ndf <- (nsmall + 2) * (nsmall - 1) / 2
1 - pchisq(mychisq, ndf)

X1 <- scale(wines[1:88, ])
X1.svd <- svd(X1)
X1.pca <- list(scores = X1.svd$u %*% diag(X1.svd$d),
               loadings = X1.svd$v)

X2 <- scale(wines[89:177, ],
            center = attr(X1, "scaled:center"),
            scale = attr(X1, "scaled:scale"))
X2.scores <- X2 %*% X1.pca$loadings

wines.odd <- seq(1, nrow(wines), by = 2)
wines.even <- seq(2, nrow(wines), by = 2)
X1b <- scale(wines[wines.odd, ])
X2b <- scale(wines[wines.even, ],
            center = attr(X1b, "scaled:center"),
            scale = attr(X1b, "scaled:scale"))

par(mfrow = c(1, 2))
labels <- rep(c(1, 2), c(89, 88))
plot(rbind(X1.pca$scores, X2.scores),
     pch = labels, col = labels,
     xlab = "PC 1", ylab = "PC 2")
X1b.svd <- svd(X1b)
X1b.pca <- list(scores = X1b.svd$u %*% diag(X1b.svd$d),
               loadings = X1b.svd$v)
X2b.scores <- X2b %*% X1b.pca$loadings
plot(rbind(X1b.pca$scores, X2b.scores),
     pch = labels, col = labels,
     xlab = "PC 1", ylab = "PC 2")

data(gasoline, package = "pls")
wavelengths <- seq(900, 1700, by = 2)
nir.prcomp <- prcomp(gasoline$NIR, rank. = 6)
summary(nir.prcomp)

par(opar)
plot(nir.prcomp, main = "Gasoline scree plot")
nir.loadings <- nir.prcomp$rotation[, 1:4]

par(opar)
par(mfrow = c(1, 2), pty = "s")
offset <- c(0, 0.09) # to create space for labels
plot(nir.loadings[, 1:2], type = "l",
     xlim = range(nir.loadings[, 1]) + offset,
     xlab = "PC 1 (72.6%)", ylab = "PC 2 (11.3%)")
points(nir.loadings[c(386, 396), 1:2])
text(nir.loadings[c(386, 396), 1:2], pos = 4,
     labels = paste(c(1670, 1690), "nm"))
offset <- c(-0.12, 0.12) # to create space for labels
plot(nir.loadings[, 3:4], type = "l",
     xlim = range(nir.loadings[, 3]) + offset,
     xlab = "PC 3 (6.95%)", ylab = "PC 4 (4.60%)")
points(nir.loadings[c(154, 370, 398), 3:4])
text(nir.loadings[c(154, 370), 3:4], pos = 4,
     labels = paste(c(1206, 1638), "nm"))
text(nir.loadings[398, 3:4, drop = FALSE],
     labels = "1694 nm", pos = 2)

par(opar)
biplot(nir.prcomp, col = c("black", "blue"), scale = 0.5)

extremes <- c(15, 41, 43, 57)
Xextr <- scale(gasoline$NIR, scale = FALSE)[extremes, ]
matplot(wavelengths, t(Xextr),
        type = "l", xlab = "Wavelength (nm)",
        ylab = "Intensity (mean-scaled)", lty = c(1, 1, 2, 2),
        col = 1:4)
legend("bottomleft", legend = paste("sample", extremes),
       lty = c(1, 1, 2, 2), col = 1:4, bty = "n")

library(MASS)
par(mfrow = c(1, 3))
wines.dist <- dist(scale(wines))
wines.cmdscale <- cmdscale(wines.dist)
plot(wines.cmdscale, 
     pch = wine.classes, col = wine.classes, 
     main = "Principal Coordinate Analysis",
     xlab = "Coord 1", ylab = "Coord 2")
wines.sammon <- sammon(wines.dist, trace = FALSE)
plot(wines.sammon$points, main = "Sammon mapping",
     col = wine.classes, pch = wine.classes, 
     xlab = "Coord 1", ylab = "Coord 2")
wines.isoMDS <- isoMDS(wines.dist, trace = FALSE)
plot(wines.isoMDS$points, main = "Non-metric MDS",
     col = wine.classes, pch = wine.classes, 
     xlab = "Coord 1", ylab = "Coord 2")

library(fastICA)
set.seed(7)
wines.ica <- fastICA(wines.sc, 3)
pairs(wines.ica$S, main = "ICA components",
      col = wine.classes, pch = wine.classes)

set.seed(7)
wines.ica5 <- fastICA(wines.sc, 5)
pairs(wines.ica5$S[, 1:3], 
      main = "ICA components (3 out of 5)",
      col = wine.classes, pch = wine.classes)

(wines.fa <- factanal(wines.sc, 3, scores = "regression"))

pairs(wines.fa$scores, main = "FA components",
      col = wine.classes, pch = wine.classes)

