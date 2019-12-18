## Demo featuring all R code from Chapter 6 in the book "Chemometrics
## with R" by Ron Wehrens. The only additions to the book chapter code
## are loading the required libraries and data, and occasionally
## closing a graphics window.
graphics.off()
opar <- par(no.readonly = TRUE)

par(mar = c(1, 1, 1, 1))
clpoints <- matrix(c(.5, .65, .85, 1,
                     4.5, 4.35, 4.15, 4,
                     .5, .65, .85, 1,
                     .5, .65, .85, 1), ncol = 2)
plot(clpoints, ylim = c(-.5, 2), xlim = c(0, 5), xlab = "", ylab = "",
     axes = FALSE, cex = 1.5, pch = rep(1:2, each = 4))
box()
points(c(.75, 4.25), c(.75, .75), pch = 4, cex = .4)
arrows(c(1, .5, .75),
       c(1, .5, .75),
       c(4, 4.5, 4.25),
       c(1, .5, .75), angle = 20, length = .15,
       code = 3, col = "gray", lty = 2)
text(rep(2.5, 3), c(1, .5, .75), adj = c(.5, 1.5),
     labels = c("Single linkage", "Complete linkage", "Average linkage"))
legend("topleft", legend = paste("Cluster", 1:2), pch = 1:2)
text(clpoints, labels = c(letters[1:4], LETTERS[1:4]),
     pos = rep(c(2, 4), each = 4))

data(wines, package = "kohonen")
wine.classes <- as.integer(vintages)
wines.sc <- scale(wines)
wines.odd <- seq(1, nrow(wines), by = 2)
wines.even <- seq(2, nrow(wines), by = 2)
twowines <- wines[vintages != "Barolo",]
twovintages <- factor(vintages[vintages != "Barolo"])

par(opar)
par(mfrow = c(2, 1))
set.seed(7)
subset <- sample(nrow(wines), 20)
wines.dist <- dist(wines.sc[subset, ])
wines.hcsingle <-  hclust(wines.dist, method = "single")
plot(wines.hcsingle, labels = vintages[subset])
wines.hccomplete <- hclust(wines.dist, method = "complete")
plot(wines.hccomplete, labels = vintages[subset])

wines.cl.single <- cutree(wines.hcsingle, h = 3)
table(wines.cl.single, vintages[subset])

wines.dist <- dist(wines.sc)
wines.hcsingle <- hclust(wines.dist, method = "single")
table(vintages, cutree(wines.hcsingle, k = 3))

wines.hccomplete <- hclust(wines.dist, method = "complete")
table(vintages, cutree(wines.hccomplete, k = 3))

library(cluster)
wines.agness <- agnes(wines.dist, method = "single")
wines.agnesa <- agnes(wines.dist, method = "average")
wines.agnesc <- agnes(wines.dist, method = "complete")

cbind(wines.agness$ac, wines.agnesa$ac, wines.agnesc$ac)

(wines.km <- kmeans(wines.sc, centers = 3))

table(vintages, wines.km$cluster)

(wines.pam <- pam(wines.dist, k = 3))

par(opar)
plot(wines.pam, main = "Silhouette plot")

best.pam <- pam(wines.dist, k = 2)
for (i in 3:10) {
  tmp.pam <- pam(wines.dist, k = i)
  if (tmp.pam$silinfo$avg.width < best.pam$silinfo$avg.width)
    best.pam <- tmp.pam
}
best.pam$medoids

table(vintages, best.pam$clustering)

require(mclust, quiet = TRUE)
wines.BIC <- mclustBIC(wines.sc, modelNames = "VVV", verbose = FALSE)
plot(wines.BIC)

wines.mclust2 <- mclustModel(wines.sc, wines.BIC)
wines.mclust3 <- mclustModel(wines.sc, wines.BIC, G = 3)

par(mfrow = c(2, 2))
coordProj(wines.sc, dimens = c(7, 13),
          parameters = wines.mclust2$parameters,
          z = wines.mclust2$z, what = "classification")
title("2 clusters: classification")
coordProj(wines.sc, dimens = c(7, 13),
          parameters = wines.mclust3$parameters,
          z = wines.mclust3$z, what = "classification")
title("3 clusters: classification")
coordProj(wines.sc, dimens = c(7, 13),
          parameters = wines.mclust2$parameters,
          z = wines.mclust2$z, what = "uncertainty")
title("2 clusters: uncertainty")
coordProj(wines.sc, dimens = c(7, 13),
          parameters = wines.mclust3$parameters,
          z = wines.mclust3$z, what = "uncertainty")
title("3 clusters: uncertainty")

wines.BIC <- mclustBIC(wines.sc, verbose = FALSE)
par(opar)
plot(wines.BIC, legendArgs = list(x = "bottomleft", ncol = 3))

AdjRkl <- function(part1, part2) {
  confusion <- table(part1, part2)

  n <- sum(confusion)
  a <- sum(choose(confusion[confusion>1], 2))
  b <- apply(confusion, 1, sum)
  b <- sum(choose(b[b>1], 2))
  c <- apply(confusion, 2, sum)
  c <- sum(choose(c[c>1], 2))

  Rexp <- b*c/choose(n, 2)
  (a - Rexp) / (.5*(b+c) - Rexp )
}

library(kohonen)
set.seed(7)
som.wines <- som(wines.sc, grid = somgrid(6, 4, "hexagonal"))
set.seed(17)
som.wines2 <- som(wines.sc, grid = somgrid(6, 4, "hexagonal"))

par(mfrow = c(1, 2))
som.hc <- cutree(hclust(dist(getCodes(som.wines, 1))), k = 3)
som.hc2 <- cutree(hclust(dist(getCodes(som.wines2, 1))), k = 3)
plot(som.wines, "mapping", bgcol = terrain.colors(3)[som.hc],
     pch = as.integer(vintages), main = "Seed 7")
plot(som.wines2, "mapping", bgcol = terrain.colors(3)[som.hc2],
     pch = as.integer(vintages), main = "Seed 17")

som.clust <- som.hc[som.wines$unit.classif]
som.clust2 <- som.hc2[som.wines2$unit.classif]
AdjRkl(som.clust, som.clust2)

AdjRkl(vintages, som.clust)
AdjRkl(vintages, som.clust2)

