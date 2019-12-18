lcmsimage <-
  function(z, rt, mz,
           xlim = range(rt),
           ylim = range(mz),
           zmin = 0,
           xlab = "Time (s)", ylab = "m/z", 
           ColourScale = c("exponential", "linear"),
           nBreaks = 12,
           colors = terrain.colors(nBreaks),
           PlotProjections = c("max", "sum", "none"),
           ncolorbarticks = 5,
           colorbarticks =
             axisTicks(zlim2, log = ColourScale == "exponential",
                       nint = ncolorbarticks),
           ...)
{
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  
  ColourScale <- match.arg(ColourScale)
  
  zorig <- z
  z <- z - zmin + 1 # eventually, the minimal value of z will be 1
  z[z < 1] <- NA
  zlim <- c(1, max(z, na.rm = TRUE))
  
  ## When there is not enough space for the color bar it is
  ## invisible. Needs at least 5 cm. Let's calculate the fraction of
  ## the total width. For the MS.wd, the same holds.
  t.wd <- par("din")[1] * 2.54
  MS.wd <- clrbar.wd <- 4 / t.wd
  
  PlotProjections <- match.arg(PlotProjections)
  
  if (PlotProjections != "none") {
    layout(matrix(c(5,1,4,3,6,2), 2, 3, byrow = TRUE), 
           heights = c(1,3),
           widths = c(clrbar.wd, 1 - clrbar.wd - MS.wd, MS.wd))
    topmar <- rightmar <- 1
  } else {
    layout(matrix(c(1,2), 1, 2, byrow = TRUE), 
           widths = c(clrbar.wd, 1 - clrbar.wd))
    topmar <- 4
    rightmar <- 2
  }
  
  ## For the marginals, we switch back to the original data
  zorig[is.na(zorig)] <- 0
  
  ## Plot the chromatogram at the top, if desired
  if (PlotProjections != "none") {
    Chrom <- apply(zorig, 1, PlotProjections)
    
    par(mar = c(1, 5, 2, 1), ann = FALSE, bty = "n", xaxs= "i")
    plot(rt, Chrom, xlim = xlim, ylim = c(0, max(Chrom)),
         type = "l", xaxs = "i", axes = FALSE, ...)
    axis(2, col = "gray", col.axis = "gray")
    abline(h = 0, col = "gray")
   
    ## Plot MS projection
    TMS <- apply(zorig, 2, PlotProjections)
    
    par(mar = c(5, 2, 1, 1), ann = FALSE, bty = "n", yaxs = "i")
    plot(TMS, mz, xlim = c(0, max(TMS)),
         type = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "",
         axes = FALSE, ...)
    abline(v = 0, col = "gray") 
    segments(0, mz, TMS, mz)
    axis(1, col = "gray", col.axis = "gray")
  }

  if (ColourScale == "exponential") {
    z <- log(z)
    z[!is.finite(z)] <- NA
    zlim2 <- log(zlim)
    if (any(!is.finite(zlim2))) {
      warning("Warning: rescaling zlim to finite values")
      if (!is.finite(zlim2[1])) zlim2[1] <- max(min(z, na.rm = TRUE), 1)
      if (!is.finite(zlim2[2])) zlim2[2] <- max(z, na.rm = TRUE)
    }
  } else {
    zlim2 <- zlim
  }

  if (length(colors) != nBreaks) nBreaks <- length(colors)
  Breaks <- seq(zlim2[1], zlim2[2], length = nBreaks+1)

  ## Draw a suitable colourbar at the left.
  par(mar = c(5, 0, topmar, 1), xaxt = "n", xpd = NA)
  plot(NA, NA, xlim = c(0, 1), ylim = zlim2, type = "n",
       ann = FALSE, bty = "n", xaxt = "n", yaxt = "n",
       xaxs = "i", yaxs = "i", ...)
  
  rect(0.9, Breaks[1:nBreaks], 1, Breaks[1+1:nBreaks], 
       col = colors, if (nBreaks > 50) {border = NA})
  rect(0.9, Breaks[1], 1, Breaks[nBreaks], lwd = 1)
  par(yaxs = "i")
  
  if (ColourScale == "exponential") {
    tm.pos <- log(colorbarticks)
    ## colorbarticks <- axTicks(2,
    ##                      axp = c(max(zlim[1] + zmin - 1, 1),
    ##                              zlim[2] + zmin - 1, 3),
    ##                      log = TRUE, nintLog = Inf)
    ##  colorbarticks <- colorbarticks[colorbarticks > zmin - 1]
    ##  tm.pos <- log(colorbarticks - zmin + 1)
    
    if (min(tm.pos) > .3) {
      tm.pos <- c(0, tm.pos)
      colorbarticks <- c(zmin, colorbarticks)
    }
  } else {
    colorbarticks <- pretty(zlim + zmin - 1)
    tm.pos <- colorbarticks - zmin + 1
  }
  
  tm.idx <- tm.pos <= zlim2[2] & tm.pos >= zlim2[1]
  text(0.82, tm.pos[tm.idx], colorbarticks[tm.idx], adj = 1)
  
  if (PlotProjections != "none") {
    par(mar = c(1, 2, 2, rightmar), ann = FALSE, bty = "n")

    plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), type = "n",
         ann = FALSE, bty = "n", xaxt = "n", yaxt = "n",
         xaxs = "i", yaxs = "i")
    text(0, 0, paste("Projection:\n", PlotProjections, sep = ""),
         adj = c(0, 0), cex = 1.2)

    plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), type = "n",
         ann = FALSE, bty = "n", xaxt = "n", yaxt = "n",
         xaxs = "i", yaxs = "i")
    text(1, 0, paste("Colour scale:\n", ColourScale, sep = ""),
         adj = c(1, 0), cex = 1.2)
  }
  
  ## Set margins for the main window
  par(mar = c(5, 5, topmar, rightmar), ann = TRUE, bty = "o")
  
  ## 2D LC-MS image
  ## Open a plot window with the correct axes
  plot(NA, NA, xlim = xlim, ylim = ylim, type = "n", xlab = xlab,
       ylab = ylab, xaxs = "i", yaxs = "i", xaxt="s", yaxt="s", ...)

  rect(xlim[1], ylim[1], xlim[2], ylim[2], col = "white")
  
  image(rt, mz, z, xlim = xlim, ylim = ylim, zlim = zlim2,
        col = colors,
        breaks = Breaks, xlab = xlab, ylab = ylab, add = TRUE, ...)
  box()
}

