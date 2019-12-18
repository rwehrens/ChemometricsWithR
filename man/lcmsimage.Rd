\name{lcmsimage}
\alias{lcmsimage}
\title{Display LC-MS data}
\description{
  LC-MS data are characterised by peak positions in two dimensions, the
  retention time dimension and the mass-to-charge ratio (m/z)
  dimension. The plot shows this plane (corresponding to one sample,
  basically an intensity matrix), with projections to the top and
  the right.
}
\usage{
lcmsimage(z, rt, mz, xlim = range(rt), ylim = range(mz),
          zmin = 0, xlab = "Time (s)", ylab = "m/z",
          ColourScale = c("exponential", "linear"), nBreaks = 12,
          colors = terrain.colors(nBreaks),
          PlotProjections = c("max", "sum", "none"),
          ncolorbarticks = 5,
          colorbarticks = axisTicks(zlim2, log = ColourScale == "exponential",
          nint = ncolorbarticks), ...)
}
\arguments{
  \item{z}{the intensity matrix}
  \item{rt}{the vector of time points corresponding to the x axis}
  \item{mz}{the vector of m/z values corresponding to the y axis}
  \item{xlim,ylim}{defining which part of the data is shown, by default all}
  \item{zmin}{minimum of the intensity range, default zero}
  \item{xlab,ylab}{axis labels}
  \item{ColourScale}{whether to use an exponential (default) or linear
    color scale} 
  \item{nBreaks,colors}{which and how many colours to use}
  \item{PlotProjections}{ways to calculate the projections to the top
    and right of the figure. Choosing \code{'none'} here suppresses the
    projections. Default is \code{'max'} (corresponding to a silhouette
    view, only the largest value is shown) but \code{'sum'} is also
    possible, which leads to the total-ion chromatogram and the direct
    injection mass spectrum.}
  \item{ncolorbarticks,colorbarticks}{control over the tick marks of the
    color bar}
  \item{\dots}{other arguments for the projection plots}
}
\author{
Ron Wehrens, Tom Bloemberg
}
\examples{
data(lcms,package = "ptw")
mycols <- terrain.colors(40)
lcmsimage(t(lcms[,,1]),
          rt = seq(2000, 5500, length = 2000),
          mz = seq(550, 599.5, length = 100),
          zmin = 1, colors = mycols, ncolorbarticks = 10)
}
\keyword{graphics}% use one of  RShowDoc("KEYWORDS")
