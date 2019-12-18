## Demo featuring all R code from Chapter 2 in the book "Chemometrics
## with R" by Ron Wehrens. The only additions to the book chapter code
## are loading the required libraries and data, and occasionally
## closing a graphics window.

data(gasoline, package = "pls")
wavelengths <- seq(900, 1700, by = 2)
matplot(wavelengths, t(gasoline$NIR), type = "l", 
        lty = 1, xlab = "Wavelength (nm)", ylab = "1/R")

data(wines, package = "kohonen")
colnames(wines)
table(vintages)

wine.classes <- as.integer(vintages)
pairs(wines[, 1:3], pch = wine.classes, col = wine.classes)

data(Prostate2000Raw)
plot(Prostate2000Raw$mz, Prostate2000Raw$intensity[, 1],
     type = "h", main = "Prostate data",
     xlab = bquote(italic(.("m/z"))~.("(Da)")),
     ylab = "Intensity")
table(Prostate2000Raw$type)

data(lcms, package = "ptw")
require(RColorBrewer)
mycols <- colorRampPalette(brewer.pal(9,"Blues"))(40)
lcmsimage(t(lcms[, , 1]),
          rt = seq(2000, 5500, length = 2000),
          mz = seq(550, 599.5, length = 100),
          zmin = 1, colors = mycols, ncolorbarticks = 10)

