## Demo featuring all R code from Chapter 10 in the book "Chemometrics
## with R" by Ron Wehrens. The only additions to the book chapter code
## are loading the required libraries and data, and occasionally
## closing a graphics window.
graphics.off()
opar <- par(no.readonly = TRUE)

require(kohonen, quiet = TRUE) # for classvec2classmat
data(wines, package = "kohonen")
wine.classes <- as.integer(vintages)
wines.sc <- scale(wines)
wines.odd <- seq(1, nrow(wines), by = 2)
wines.even <- seq(2, nrow(wines), by = 2)
twowines <- wines[vintages != "Barolo",]
twovintages <- factor(vintages[vintages != "Barolo"])

X <- wines[wines.odd, ]                     
C <- classvec2classmat(vintages[wines.odd])
wines.lm <- lm(C ~ X)                
wines.lm.summ <- summary(wines.lm)                    
wines.lm.summ[[3]]

sapply(wines.lm.summ, 
       function(x) which(x$coefficients[, 4] < .1))

require(boot, quiet = TRUE)
require(pls, quiet = TRUE)
data(gasoline, package = "pls")
wavelengths <- seq(900, 1700, by = 2)
gas.odd <- seq(1, nrow(gasoline$NIR), by = 2)
gas.even <- seq(2, nrow(gasoline$NIR), by = 2)

## The following object is also used in chapter9.R
npc <- 4
set.seed(17)
gas.pcr.bootCI <-
  boot(gasoline,
       function(x, ind) {
         c(coef(pcr(octane ~ ., data = x,
                    ncomp = npc, subset = ind)))},
       R = 999)
set.seed(7)
gas.BCACI <-
  t(sapply(1:ncol(gasoline$NIR),
           function(i, x) {
             boot.ci(x, index = i, type = "bca")$bca[, 4:5]},
           gas.pcr.bootCI))

BCAcoef <- gas.pcr.bootCI$t0
signif <- gas.BCACI[, 1] > 0 | gas.BCACI[, 2] < 0
BCAcoef[!signif] <- NA
matplot(wavelengths, gas.BCACI, type = "n",
        xlab = "Wavelength (nm)",
        ylab = "Regression coefficient",
        main = "Gasoline data: PCR (4 PCs)")
abline(h = 0, col = "gray")
polygon(c(wavelengths, rev(wavelengths)),
        c(gas.BCACI[, 1], rev(gas.BCACI[, 2])),
        col = "pink", border = NA)
lines(wavelengths, BCAcoef, lwd = 2)

smallmod <- pcr(octane ~ NIR[, signif], data = gasoline,
                ncomp = 4, validation = "LOO")
RMSEP(smallmod, intercept = FALSE, estimate = "CV")

twowines.df <- data.frame(vintage = twovintages, twowines)
twowines.lm0 <- lm(as.integer(vintage) ~ 1, data = twowines.df)
add1(twowines.lm0, scope = names(twowines.df)[-1])

twowines.lmfull <- lm(as.integer(vintage) ~ ., data = twowines.df)
drop1(twowines.lmfull)

step(twowines.lmfull, trace = 0)

require(leaps, quiet = TRUE)
twowines.leaps <- regsubsets(vintage ~ ., data = twowines.df)
twowines.leaps.sum <- summary(twowines.leaps)
names(which(twowines.leaps.sum$which[8, ]))

twowines <- twowines.df[, -1]
vint <- twowines.df[, 1]

wines.counts <- table(vint)
wines.groups <- split(as.data.frame(twowines), vint)
wines.covmats <- lapply(wines.groups, cov)
wines.wcovmats <- lapply(1:length(wines.groups),
                         function(i, x, y) x[[i]]*y[i],
                         wines.covmats, wines.counts)
wines.pcov12 <- Reduce("+", wines.wcovmats) / (length(vint) - 2)
MLLDA <-
  solve(wines.pcov12,
        apply(sapply(wines.groups, colMeans), 1, diff))
WSS <- Reduce("+", lapply(wines.groups,
                          function(x) {
                            crossprod(scale(x, scale = FALSE))}))
BSS <- Reduce("+", lapply(wines.groups,
                          function(x, y) {
                            nrow(x) * tcrossprod(colMeans(x) - y)},
                          colMeans(twowines)))

Tii <- solve(BSS + WSS)
Ddist <- mahalanobis(colMeans(wines.groups[[1]]),
                     colMeans(wines.groups[[2]]),
                     wines.pcov12)
m <- sum(sapply(wines.groups, nrow)) - 2
p <- ncol(wines)
c <- prod(sapply(wines.groups, nrow)) / 
  sum(sapply(wines.groups, nrow))
Fcal <- (MLLDA^2 / diag(Tii)) * 
  (m - p + 1) * c^2 / (m * (m + c^2 * Ddist))
which(Fcal > qf(.95, 1, m-p+1))

require(glmnet, quiet = TRUE)
par(mfrow = c(1, 2))
gas.lasso <- glmnet(x = gasoline$NIR[gas.odd, ],
                    y = gasoline$octane[gas.odd])
plot(gas.lasso, xvar = "lambda", label = TRUE)
set.seed(7)
gas.lasso.cv <- cv.glmnet(gasoline$NIR[gas.odd, ], 
                          gasoline$octane[gas.odd])
svals <- gas.lasso.cv[c("lambda.1se", "lambda.min")]
abline(v = svals, col = "darkgray", lty = 2)
plot(gas.lasso.cv)

gas.lasso.preds <-
  lapply(svals,
         function(x)
           predict(gas.lasso,
                   newx = gasoline$NIR[gas.even, ],
                   s = x))
sapply(gas.lasso.preds,
       function(x) rms(x, gasoline$octane[gas.even]))

gas.lasso.coefs <- lapply(svals,
                          function(x) coef(gas.lasso, s = x))
sapply(gas.lasso.coefs,
       function(x) sum(x != 0))

set.seed(7)
gas.elnet.cv <- cv.glmnet(gasoline$NIR, gasoline$octane, alpha = .5)
elnet.svals <- c("one.se" = gas.elnet.cv$lambda.1se,
                 "min" = gas.elnet.cv$lambda.min)
snames <- c("lambda.1se", "lambda.min") 
gas.elnet.preds <-
  lapply(snames,
         function(x)
           predict(gas.elnet.cv,
                   newx = gasoline$NIR[gas.even, ],
                   s = x))
names(gas.elnet.preds) <- snames

gas.elnet.coefs <- lapply(snames,
                          function(x) coef(gas.elnet.cv, s = x))

par(mfrow = c(1, 2))
gas.elnet <- glmnet(gasoline$NIR, gasoline$octane, alpha = .5)
plot(gas.elnet, "norm")
set.seed(7)
gas.elnet.cv <- cv.glmnet(gasoline$NIR, gasoline$octane, alpha = .5)
elnet.svals <- c("one.se" = gas.elnet.cv$lambda.1se,
                 "min" = gas.elnet.cv$lambda.min)
snames <- c("lambda.1se", "lambda.min") 
gas.elnet.preds <-
  lapply(snames,
         function(x)
           predict(gas.elnet.cv,
                   newx = gasoline$NIR[gas.even, ],
                   s = x))
names(gas.elnet.preds) <- snames

gas.elnet.coefs <- lapply(snames,
                          function(x) coef(gas.elnet.cv, s = x))
abline(v = elnet.svals, col = "darkgray", lty = 2)
plot(gas.elnet.cv)

sapply(gas.elnet.preds,
       function(x) rms(x, gasoline$octane[gas.even]))
sapply(gas.elnet.coefs,
       function(x) sum(x != 0))

par(opar)
offset <- 15
plot(wavelengths, gas.lasso.coefs[[2]][-1], type = "h",
     ylab = "Coefficients", xlab = "Wavelength (nm)", lwd = 3)
huhn <- par("usr")
par(usr = huhn + c(0, 0, -offset, -offset))
lines(wavelengths, gas.elnet.coefs[[2]][-1], type = "h",
      lwd = 2, col = "blue")
axis(4, col = "blue", col.axis = "blue")
legend("topleft", lty = 1, col = c("black", "blue"),
       legend = c("Lasso", "Elastic net"), bty = "n")

require(MASS, quiet = TRUE)
lda.loofun <- function(selection, xmat, grouping, ...) {
  if (sum(selection) == 0) return(100)
  lda.obj <- lda(xmat[, selection == 1], grouping, CV = TRUE)
  100*sum(lda.obj$class != grouping)/length(grouping)
}

saStepFun <- function(selected, ...) {
  maxval <- length(selected)
  selection <- which(selected == 1)
  newvar2 <- sample(1:maxval, 2)

  ## too short: add a random number
  if (length(selection) < 2) {
    result <- unique(c(selection, newvar2))[1:2]
  } else {     # generate two variable numbers 
    presentp <- newvar2 %in% selection
    ## if both are in x, remove the first
    if (all(presentp)) {
      result <- selection[selection != newvar2[1]]
    } else {   # if none are in selection, add the first
      if (all(!presentp)) {
        result <- c(selection, newvar2[1])
      } else { # otherwise swap
        result <- c(selection[selection != newvar2[presentp]],
                    newvar2[!presentp])
      }}}
  
  newselected <- rep(0, length(selected))
  newselected[result] <- 1
  newselected
}

set.seed(173)
initselect <- rep(0, ncol(wines))
initselect[sample(1:ncol(wines), 5)] <- 1
(r0 <- lda.loofun(initselect, x = twowines,
                  grouping = twovintages))

set.seed(7)
SAoptimWines <-
  optim(initselect,
        fn = lda.loofun, gr = saStepFun, method = "SANN", 
        x = twowines, grouping = twovintages)

SAoptimWines[c("par", "value")]

pls.cvfun <- function(selection, xmat, response, ncomp, ...) {
  if (sum(selection) < ncomp) return(Inf)
  pls.obj <- plsr(response ~ xmat[, selection == 1],
                  validation = "CV", ncomp = ncomp, ...)
  c(RMSEP(pls.obj, estimate = "CV", ncomp = ncomp,
         intercept = FALSE)$val)
}

set.seed(7)
nNIR <- ncol(gasoline$NIR)
initselect <- rep(0, nNIR)
initselect[sample(1:nNIR, 8)] <- 1
SAoptimNIR1 <-
  optim(initselect,
        fn = pls.cvfun, gr = saStepFun, method = "SANN",
        x = gasoline$NIR, response = gasoline$octane,
        ncomp = 2, maxval = nNIR)

pls.cvfun(initselect, gasoline$NIR, gasoline$octane, ncomp = 2)
(nvarSA1 <- sum(SAoptimNIR1$par))
SAoptimNIR1$value

pls.cvfun2 <- function(selection, xmat, response, ncomp,
                       penalty = 0.01, ...) {
  if (sum(selection) < ncomp) return(Inf)
  pls.obj <- plsr(response ~ xmat[, selection == 1, drop = FALSE],
                  validation = "CV", ncomp = ncomp, ...)
  c(RMSEP(pls.obj, estimate = "CV", ncomp = ncomp,
         intercept = FALSE)$val) +
    penalty * sum(selection)
}

saStepFun2 <- function(selected, plimits = c(.3, .7), ...) {
  dowhat <- runif(1)
  
  ## decrease selection
  if (dowhat < plimits[1]) { 
    if (sum(selected) > 2) { # not too small...
      kickone <- sample(which(selected == 1), 1)
      selected[kickone] <- 0
      return(selected)
    }
  }
  
  ## increase selection
  if (dowhat > plimits[2]) { # not too big...
    if (sum(selected) < length(selected)) {
      addone <- sample(which(selected == 0), 1)
      selected[addone] <- 1
      return(selected)
    }  
  }
  
  ## swap
  kickone <- sample(which(selected == 1), 1)
  selected[kickone] <- 0
  addone <- sample(which(selected == 0), 1)
  selected[addone] <- 1
  selected
}

set.seed(7)
penalty <- 0.01
SAoptimNIR2 <-
  optim(initselect,
        fn = pls.cvfun2, gr = saStepFun2, method = "SANN",
        x = gasoline$NIR, response = gasoline$octane,
        ncomp = 2, maxval = nNIR,
        control = list(temp = 1))

(nvarSA2 <- sum(SAoptimNIR2$par))
SAoptimNIR2$value - penalty*nvarSA2

require(subselect, quiet = TRUE)
set.seed(7)
winesHmat <- ldaHmat(twowines.df[, -1], twowines.df[, 1])
wines.anneal <- 
  anneal(winesHmat$mat, kmin = 3, kmax = 3,
         H = winesHmat$H, criterion = "ccr12", r = 1)

wines.anneal$bestsets
wines.anneal$bestvalues

ccr12.coef((nrow(twowines.df) - 1) * var(twowines.df[, -1]),
           winesHmat$H, r = 1, c(7, 10, 11))

selection <- rep(0, ncol(twowines))
selection[c(2, 7, 10)] <- 1
lda.loofun(selection, twowines.df[, -1], twowines.df[, 1])

fitnessfun <- function(...) -pls.cvfun2(...)

require(GA, quiet = TRUE)
gaControl("binary" = list(mutation = "gabin_raMutation"))
set.seed(7)
GAoptimNIR1 <-
  ga(type = "binary", fitness = fitnessfun,
     x = gasoline$NIR, response = gasoline$octane,
     ncomp = 2, penalty = penalty,
     nBits = ncol(gasoline$NIR), monitor = FALSE, maxiter = 100)
myMutate <- function (object, parent, bias = 0.01) 
{
  mutate <- parent <- as.vector(object@population[parent, ])
  n <- length(parent)
  probs <- abs(mutate - bias)
  j <- sample(1:n, size = 1, prob = probs)
  mutate[j] <- abs(mutate[j] - 1)
  mutate
}
set.seed(7)
gaControl("binary" = list(mutation = "myMutate"))
popSize <- 50 # default
initmat <- matrix(0, popSize, nNIR)
initmat[sample(1:(popSize*nNIR), nNIR)] <- 1

GAoptimNIR2 <-
  ga(type = "binary", fitness = fitnessfun,
     x = gasoline$NIR, response = gasoline$octane,
     popSize = popSize, nBits = ncol(gasoline$NIR),
     ncomp = 2, suggestions = initmat, penalty = penalty,
     monitor = FALSE, maxiter = 100)

(nvarGA1 <- sum(GAoptimNIR1@solution))
-GAoptimNIR1@fitnessValue + penalty*nvarGA1

(nvarGA2 <- sum(GAoptimNIR2@solution))
-GAoptimNIR2@fitnessValue + penalty*nvarGA2

par(mfrow = c(1, 2))
plot(GAoptimNIR1, main = "Naive")
plot(GAoptimNIR2, main = "Dedicated")

wines.genetic <-
  genetic(winesHmat$mat, kmin = 3, kmax = 5, nger = 20,
          popsize = 50, maxclone = 0,
          H = winesHmat$H, criterion = "ccr12", r = 1)
wines.genetic$bestvalues
wines.genetic$bestsets

