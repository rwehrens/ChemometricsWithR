## Demo featuring all R code from Chapter 7 in the book "Chemometrics
## with R" by Ron Wehrens. The only additions to the book chapter code
## are loading the required libraries and data, and occasionally
## closing a graphics window.
graphics.off()
opar <- par(no.readonly = TRUE)

data(wines, package = "kohonen")
wine.classes <- as.integer(vintages)
wines.sc <- scale(wines)
wines.odd <- seq(1, nrow(wines), by = 2)
wines.even <- seq(2, nrow(wines), by = 2)
twowines <- wines[vintages != "Barolo",]
twovintages <- factor(vintages[vintages != "Barolo"])

wines.trn <- wines[wines.odd, ]
wines.tst <- wines[wines.even, ]
vint.trn <- vintages[wines.odd]
vint.tst <- vintages[wines.even]

wines2.trn <- wines.trn[, c(7, 13)]
wines2.tst <- wines.tst[, c(7, 13)]

wines.counts <- table(vint.trn)
ngroups <- length(wines.counts)
wines.groups <- split(as.data.frame(wines.trn), vint.trn)
wines.covmats <- lapply(wines.groups, cov)
wines.wcovmats <- mapply('*', wines.covmats, wines.counts,
                         SIMPLIFY = FALSE)
wines.pooledcov <-
  Reduce("+", wines.wcovmats) / (nrow(wines.trn) - ngroups)

wines.pooledcov2 <- matrix(0, ncol(wines), ncol(wines))
for (i in 1:3) {
  wines.pooledcov2 <- wines.pooledcov2 +
    cov(wines.groups[[i]]) * nrow(wines.groups[[i]])
}
wines.pooledcov2 <- 
  wines.pooledcov2 / (nrow(wines.trn) - ngroups)

range(wines.pooledcov2 - wines.pooledcov)

distances <- 
  sapply(1:ngroups, 
         function(i, samples, means, covs)
           mahalanobis(samples, colMeans(means[[i]]), covs),
         wines.trn, wines.groups, wines.pooledcov)
trn.pred <- apply(distances, 1, which.min)

table(vint.trn, trn.pred)

distances <- 
  sapply(1:ngroups,
         function(i, samples, means, covs)
           mahalanobis(samples, colMeans(means[[i]]), covs),
         wines.tst, wines.groups, wines.pooledcov)
tst.pred <- apply(distances, 1, which.min)
table(vint.tst, tst.pred)

library(MASS)
wines.ldamod <- lda(wines.trn, grouping = vint.trn,
                    prior = rep(1, 3)/3)
wines.lda.testpred <- predict(wines.ldamod, new = wines.tst)
table(vint.tst, wines.lda.testpred$class)

par(opar)
plot(wines.ldamod, col = as.integer(vint.trn))

wines.ldamod <- lda(wines.trn, grouping = vint.trn,
                    prior = rep(1, 3)/3, CV = TRUE)
table(vint.trn, wines.ldamod$class)

wns <- wines[vintages != "Barolo", c(7, 13)]
vnt <- factor(vintages[vintages != "Barolo"])
wines.odd2 <- seq(1, nrow(wns), by = 2)
wines.even2 <- seq(2, nrow(wns), by = 2)

wines.counts <- table(vnt)
wines.groups <- split(as.data.frame(wns), vnt)
WSS <- 
  Reduce("+", lapply(wines.groups,
                     function(x) {
                       crossprod(scale(x, scale = FALSE))}))
BSS <- 
  Reduce("+", lapply(wines.groups,
                     function(x, y) {
                       nrow(x) * tcrossprod(colMeans(x) - y)},
                     colMeans(wns)))
FLDA <- eigen(solve(WSS, BSS))$vectors[, 1]
FLDA / FLDA[1]

wines.covmats <- lapply(wines.groups, cov)
wines.wcovmats <- lapply(1:length(wines.groups),
                         function(i, x, y) x[[i]]*y[i],
                         wines.covmats, wines.counts)
wines.pcov12 <- Reduce("+", wines.wcovmats) / (length(vnt) - 2)
MLLDA <-
  solve(wines.pcov12,
        apply(sapply(wines.groups, colMeans), 1, diff))
MLLDA / MLLDA[1]

wns <- wines2.trn
vnt <- vint.trn

wines.counts <- table(vnt)
wines.groups <- split(as.data.frame(wns), vnt)
WSS <- 
  Reduce("+", lapply(wines.groups,
                     function(x) {
                       crossprod(scale(x, scale = FALSE))}))
BSS <- 
  Reduce("+", lapply(wines.groups,
                     function(x, y) {
                       nrow(x) * tcrossprod(colMeans(x) - y)},
                     colMeans(wns)))
FLDA <- eigen(solve(WSS, BSS))$vectors[, 1]
FLDA / FLDA[1]

xcoo <- seq(.4, 5.4, length = 251)
ycoo <- seq(250, 1750, length = 251)
gridXY <- data.matrix(expand.grid(xcoo, ycoo))
scores <- gridXY %*% FLDA
meanscores <- c(t(sapply(wines.groups, colMeans)) %*% FLDA)
Fdistance <- outer(scores, meanscores, 
                   FUN = function(x, y) abs(x - y))
Fclassif <- apply(Fdistance, 1, which.min)

wns <- wines[vintages != "Barolo", c(7, 13)]
vnt <- factor(vintages[vintages != "Barolo"])
wines.odd2 <- seq(1, nrow(wns), by = 2)
wines.even2 <- seq(2, nrow(wns), by = 2)

par(mfrow = c(1, 2))
softbrg <- colorRampPalette(c("lightgray", "pink", "lightgreen"))
image(x = xcoo, y = ycoo,
      z = matrix(Fclassif, nrow = length(xcoo)),
      xlab = "flavonoids", ylab = "proline",
      main = "Fisher LDA", col = softbrg(3))
box()
contour(x = xcoo, y = ycoo, drawlabels = FALSE,
        z = matrix(Fclassif, nrow = length(xcoo)),
        add = TRUE)
points(wines2.tst, col = wine.classes[wines.even],
       pch = wine.classes[wines.even])
wines.ldamod <- lda(wines2.trn, 
                    grouping = vint.trn,
                    prior = rep(1, 3)/3)
colnames(gridXY) <- colnames(wines)[c(7, 13)]
lda.2Dclassif <- predict(wines.ldamod, newdata = gridXY)$class
lda.2DCM <- matrix(as.integer(lda.2Dclassif), nrow = length(xcoo))
image(x = xcoo, y = ycoo, z = lda.2DCM,
      xlab = "flavonoids", ylab = "proline",
      main = "LDA", col = softbrg(3))
box()
contour(x = xcoo, y = ycoo, z = lda.2DCM, drawlabels = FALSE,
        add = TRUE)
points(wines2.tst, col = wine.classes[wines.even],
       pch = wine.classes[wines.even])

wines2.groups <- split(as.data.frame(wines2.trn), vint.trn)
wines2.covmats <- lapply(wines2.groups, cov)
ngroups <- length(wines2.groups)
distances <- sapply(1:ngroups,
                    function(i, samples, means, covs) {
                      mahalanobis(samples,
                                  colMeans(means[[i]]),
                                  covs[[i]]) },
                    wines2.tst, wines2.groups, wines2.covmats)
test.pred <- apply(distances, 1, which.min)
table(vint.tst, test.pred)

par(mfrow = c(1, 2))
qda.mahal.dists <-
  sapply(1:ngroups,
         function(i, samples, means, covs) {
           mahalanobis(samples,
                       colMeans(means[[i]]),
                       covs[[i]]) },
         gridXY, wines2.groups, wines2.covmats)
qda.2Dclassif <- apply(qda.mahal.dists, 1, which.min)
qda.2DCM <- matrix(qda.2Dclassif, nrow = length(xcoo))
image(x = xcoo, y = ycoo, z = qda.2DCM,
      xlab = "flavonoids", ylab = "proline",
      main = "QDA", col = softbrg(3))
box()
contour(x = xcoo, y = ycoo, z = qda.2DCM, drawlabels = FALSE,
        add = TRUE)
points(wines2.tst, col = as.integer(vint.tst),
       pch = as.integer(vint.tst))

wines.qda <- qda(wines.trn, vint.trn, 
                 prior = rep(1, 3)/3)
test.qdapred <- predict(wines.qda, newdata = wines.tst)
table(vint.tst, test.qdapred$class) 

require(mclust, quiet = TRUE)
wines.MclustDA <-
  MclustDA(wines.trn, vint.trn, G = 1:5, verbose = FALSE)

summary(wines.MclustDA)

wines.McDApred <-
  predict(wines.MclustDA, newdata = wines.tst)$classification
sum(wines.McDApred != vint.tst)

wines.mclust2D <-
  MclustDA(wines2.trn, vint.trn, G = 1:5, verbose = FALSE)
wines.mclust2Dpred <- predict(wines.mclust2D, gridXY)
mbda.2DCM <- matrix(as.integer(wines.mclust2Dpred$classification),
                    nrow = length(xcoo))
image(x = xcoo, y = ycoo, z = mbda.2DCM,
      main = "MBDA", xlab = "flavonoids", ylab = "proline",
      col = softbrg(3))
box()
contour(x = xcoo, y = ycoo, z = mbda.2DCM, drawlabels = FALSE,
        add = TRUE)
points(wines2.tst,
       col = as.integer(vint.tst),
       pch = as.integer(vint.tst))

library(sfsmisc) # for dDA
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

prost.dlda <-
  dDA(prost[prost.odd, ], as.integer(prost.type)[prost.odd])

prost.dldapred <- predict(prost.dlda, prost[prost.even, ])
table(prost.type[prost.even], prost.dldapred)
rslt <- table(prost.type[prost.even], prost.dldapred)
rsltPerc <- round(100*sum(diag(rslt)) / sum(rslt))

prost.dqda <-
  dDA(prost[prost.odd, ], as.integer(prost.type)[prost.odd],
      pool = FALSE)
prost.dqdapred <- predict(prost.dqda, prost[prost.even, ])
table(prost.type[prost.even], prost.dqdapred)

library(rda)
prost.rda <-
  rda(t(prost[prost.odd, ]), as.integer(prost.type)[prost.odd],
      delta = seq(0, .4, len = 5), alpha = seq(0, .4, len = 5))

prost.rda

prost.rdacv <-                                                               
  rda.cv(prost.rda, t(prost[prost.odd, ]),
         as.integer(prost.type)[prost.odd])

prost.rdapred <-
  predict(prost.rda,
          t(prost[prost.odd, ]), as.integer(prost.type)[prost.odd],
          t(prost[prost.even, ]), alpha = .2, delta = 0)
table(prost.type[prost.even], prost.rdapred)

wines.trn.sc <- scale(wines.trn)
wines.tst.sc <- scale(wines.tst,
                      scale = apply(wines.trn, 2, sd), 
                      center = colMeans(wines.trn))

dist2sample2a <- mahalanobis(wines.trn.sc,
                             wines.tst.sc[1, ], 
                             diag(13))
dist2sample2b <- mahalanobis(wines.trn,
                             wines.tst[1, ],
                             diag(apply(wines.trn, 2, var)))

range(dist2sample2a - dist2sample2b)

nearest.classes <- vint.trn[order(dist2sample2a)]
table(nearest.classes[1:10])

dist2sample2 <- mahalanobis(wines.trn,
                            wines.tst[1, ],
                            cov(wines.trn))
nearest.classes <- vint.trn[order(dist2sample2)]
table(nearest.classes[1:10])

library(class)
set.seed(7)
knn(wines.sc[wines.odd, ], wines.sc[68, ], cl = vint.trn, k = 4)
knn(wines.sc[wines.odd, ], wines.sc[68, ], cl = vint.trn, k = 4)

knn(wines.sc[wines.odd, ], wines.sc[68, ], cl = vint.trn,
    k = 4, l = 3)

set.seed(7)
wines.knnresult <- rep(0, 10)
for (i in 1:10) {
  wines.knncv <- knn.cv(wines.sc[wines.odd, ], vint.trn, k = i)
  wines.knnresult[i] <-
    sum(diag(table(vint.trn, wines.knncv)))
}
round(100 * wines.knnresult / length(wines.odd), 1)

library(e1071)
set.seed(7)
(knn.tuned <- tune.knn(wines.sc[wines.odd, ], vint.trn, k = 1:10))

par(mfrow = c(1, 2))
plot(knn.tuned)

set.seed(7)
bestKs <- rep(0, 1000)
for (i in 1:1000)
  bestKs[i] <- tune.knn(wines.sc[wines.odd, ],
                        vint.trn, 
                        k = 1:10)$best.parameters[1, 1]
hist(bestKs, breaks = 0:10 + .5)

library(rpart)
wines2.df <- data.frame(vint = vintages, wines[, c(7, 13)])
wines2.rpart <- rpart(vint ~ ., subset = wines.odd,
                      data = wines2.df, method = "class")

wines2.rpart

par(mfrow = c(1, 2))
plot(wines2.rpart, margin = .12)
text(wines2.rpart, use.n = TRUE)
cl.id <- as.integer(wines2.df$vint[wines.odd])
plot(wines2.df[wines.odd, 2:3], pch = cl.id, col = cl.id)
segments(wines2.rpart$splits[1, 4], par("usr")[3],
         wines2.rpart$splits[1, 4], par("usr")[4], lty = 2)
segments(wines2.rpart$splits[1, 4], wines2.rpart$splits[2, 4],
         par("usr")[2], wines2.rpart$splits[2, 4], lty = 2)

wines.df <- data.frame(vint = vintages, wines)
wines.rpart <- rpart(vint ~ ., subset = wines.odd,
                     data = wines.df, method = "class")

par(opar)
par(mar = c(0, 0, 0, 0))
plot(wines.rpart, margin = 0.15, compress = TRUE)
text(wines.rpart, use.n = TRUE)

wines.rpart.predict <- predict(wines.rpart, 
                               newdata = wines.df[wines.even, ])
wines.rpart.predict[31:34, ]

par(opar)
matplot(wines.rpart.predict, xlab = "sample number (test set)")

table(vint.tst, 
      predict(wines.rpart, newdata = wines.df[wines.even, ], 
              type = "class"))

gini <- function(x, clss) {
  p <- table(clss) / length(clss)
  gini.parent <- 1 - sum(p^2)
  
  gini.index <- 
    sapply(sort(x), function(splitpoint) {
      left.ones <- clss[x < splitpoint]  
      right.ones <- clss[x >= splitpoint]
      nleft <- length(left.ones)                                        
      nright <- length(right.ones)
      
      if ((nleft == 0) | (nright == 0)) return (NA)

      p.left <- table(left.ones) / nleft
      p.right <- table(right.ones) / nright

      (nleft * (1 - sum(p.left^2)) +
       nright * (1 - sum(p.right^2))) /
        (nleft + nright)
    })
  
  gini.index - gini.parent
}

wines2.df.odd <- wines2.df[wines.odd, ] 
Ginis <- sapply(wines2.df.odd[, -1], gini, wines2.df.odd$vint)
apply(Ginis, 2, min, na.rm = TRUE)
(idx <- which.min(Ginis[, 1]))
(bestSplit <- sort(wines2.df.odd[, "flavonoids"])[idx])

wl <- wines2.df.odd[wines2.df.odd$flavonoids < bestSplit, ]
GinisL <- sapply(wl[, -1], gini, wl$vint)
apply(GinisL, 2, min, na.rm = TRUE)
(idxL <- which.min(GinisL[, 2]))
(bestSplitL <- sort(wl[, "proline"])[idxL])

wr <- wines2.df.odd[wines2.df.odd$flavonoids >= bestSplit, ]
GinisR <- sapply(wr[, -1], gini, wr$vint)
apply(GinisR, 2, min, na.rm = TRUE)
(idxR <- which.min(GinisR[, 2]))
(bestSplitR <- sort(wr[, "proline"])[idxR])

par(mfrow = c(1, 2))
matplot(Ginis, pch = 1:2, col = 1:2)
legend("top", legend = colnames(Ginis),
       col = 1:2, pch = 1:2)
matplot(GinisR, pch = 1:2, col = 1:2)

prost.df <- data.frame(type = prost.type, prost = prost)
prost.rpart <- 
  rpart(type ~ ., data = prost.df, subset = prost.odd,
        control = rpart.control(cp = 0, minsplit = 0))

prost.rpartpred <- 
  predict(prost.rpart, newdata = prost.df[prost.even, ])
table(prost.type[prost.even], classmat2classvec(prost.rpartpred))

printcp(prost.rpart)

par(opar)
plotcp(prost.rpart)

prost.rpart2 <- 
  rpart(type ~ ., data = prost.df, subset = prost.odd,
        control = rpart.control(cp = 0.12))
prost.rpart2pred <- 
  predict(prost.rpart2, newdata = prost.df[prost.even, ])
table(prost.type[prost.even], classmat2classvec(prost.rpart2pred))

set.seed(7)
Sigma <- matrix(c(6, 4, 4, 3), 2, 2)/5
x1 <- mvrnorm(10, c(0, 0), Sigma, empirical = TRUE)
x2 <- mvrnorm(10, c(5, 0), Sigma, empirical = TRUE)
plot(rbind(x1, x2),
     col = c(rep(2, 10), rep(4, 10)),
     pch = c(rep(2, 10), rep(3, 10)),
     ann = FALSE)

### Horizontal class boundary
abline(v = c(max(x1[, 1]),
         mean(c(min(x2[, 1]), max(x1[, 1]))),
         min(x2[, 1])),
       lty = c(2, 1, 2))
### Class boundary with slope 2/3
margin2 <- max(x2[, 2] - x2[, 1] * 2/3)
margin1 <- min(x1[, 2] - x1[, 1] * 2/3)
abline(margin2, 2/3, lty = 2)
abline(margin1, 2/3, lty = 2)
abline(mean(c(margin1, margin2)), 2/3)

wns.df <- 
  data.frame(vint = vnt,
             flavonoids = wns[,"flavonoids"], 
             proline = wns[,"proline"])
wns.svm <- svm(vint ~ ., data = wns.df[wines.odd2, ])
wns.svmpred <- predict(wns.svm, wns.df[wines.even2, ])
table(wns.df$vint[wines.even2], wns.svmpred)

set.seed(7)
prost.svm <- svm(type ~ ., data = prost.df, subset = prost.odd,
                 cross = 10)
summary(prost.svm)

prost.svmpred <- predict(prost.svm, newdata = prost.df[prost.even, ])
table(prost.type[prost.even], prost.svmpred)

wines.svm <- svm(vint ~ flavonoids + proline, data = wines.df,
                 subset = wines.odd)
wines.svmpred.trn <- predict(wines.svm)
wines.svmpred.tst <-
  predict(wines.svm, newdata = wines.df[wines.even, ])
sum(wines.svmpred.trn == vint.trn) / length(wines.odd)
sum(wines.svmpred.tst == vint.tst) / length(wines.even)

set.seed(7)
wines.bestsvm <-
  best.svm(vint ~ flavonoids + proline, data = wines.df,
           kernel = "polynomial",
           coef0 = seq(-.5, .5, by = .1),
           gamma = 2^(-1:1), cost = 2^(2:4))

wines.bestsvmpred.trn <- 
  predict(wines.bestsvm, newdata = wines.df[wines.odd, ])
wines.bestsvmpred.tst <- 
  predict(wines.bestsvm, newdata = wines.df[wines.even, ])
sum(wines.bestsvmpred.trn == vint.trn) / length(vint.trn)
sum(wines.bestsvmpred.tst == vint.tst) / length(vint.tst)

plot(wines.svm, wines.df[wines.even, ], proline ~ flavonoids,
     color.palette = softbrg)
plot(wines.bestsvm, data = wines.df[wines.even, ],
     proline ~ flavonoids, color.palette = softbrg)

par(opar)
par(mar = rep(.2, 4))
xvals <- c(rep(1, 3), 1.5, rep(3, 4), 3.5, rep(5, 2))
yvals <- c(1:3, 4, 1:4 - .5, 4, 1:2 + .5)
plot(xvals, yvals, type = "n", axes = FALSE, ann = FALSE,
     xlim = c(0.4, 5.6), ylim = c(0.4, 4.4))
symbols(xvals, yvals,
        circles = rep(.25, length(xvals)),
        inches = FALSE, add = TRUE)
## input to hidden
arrows(rep(xvals[1:3], each = 4) + .25,
       rep(yvals[1:3], each = 4),
       rep(xvals[5:8], 3) - .25,
       rep(yvals[5:8], 3),
       angle = 15, length = .1)
## hidden to out
arrows(rep(xvals[5:8], each = 2) + .25,
       rep(yvals[5:8], each = 2),
       rep(xvals[10:11], 4) - .25,
       rep(yvals[10:11], 4),
       angle = 15, length = .1)
## bias units
arrows(rep(xvals[4], each = 4) + .25,
       rep(yvals[4], each = 4),
       xvals[5:8] - .25,
       yvals[5:8],
       angle = 15, length = .1)
arrows(rep(xvals[9], each = 2) + .25,
       rep(yvals[9], each = 2),
       xvals[10:11] - .25,
       yvals[10:11],
       angle = 15, length = .1)
text(xvals[c(4, 9)], yvals[c(4, 9)], c("bias", "bias"),
     adj = c(.5, .5))

membership.trn <- classvec2classmat(vint.trn)
head(membership.trn, 3)

library(nnet)
set.seed(7)
wines.nnet <- nnet(x = wines.trn.sc,
                   y = membership.trn,
                   size = 4)

membership.pred <- predict(wines.nnet)
training.pred <- classmat2classvec(membership.pred)
table(vint.trn, training.pred)

table(vint.tst,
      classmat2classvec(predict(wines.nnet, wines.tst.sc)))

set.seed(7)
wines.trn.sc.df <- data.frame(vintage = vint.trn, wines.trn.sc)
(wines.nnetmodels <- 
   tune.nnet(vintage ~ ., data = wines.trn.sc.df,
             size = 1:8, trace = FALSE))

par(opar)
plot(wines.nnetmodels)

set.seed(7)
best.wines.nnet <-
  best.nnet(vintage ~ ., data = wines.trn.sc.df,
            size = 1:8, trace = FALSE)
table(vint.tst,
      predict(best.wines.nnet,
              newdata = data.frame(wines.tst.sc),
              type = "class"))

