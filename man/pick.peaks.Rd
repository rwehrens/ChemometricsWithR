\name{pick.peaks}
\alias{pick.peaks}
\title{Peak-picking function.}
\description{
Function to identify local maxima in a vector, typically a spectrum or a
chromatogram.}
\usage{
pick.peaks(x, span)
}
\arguments{
  \item{x}{Numerical vector.}
  \item{span}{Neighbourhood, used to define local maxima.}
}
\value{A vector containing positions of local maxima in the input data.}
\author{Ron Wehrens}
\examples{
if (require("ptw")) {
data(lcms, package = "ptw")
plot(lcms[1,,1], type = "l", xlim = c(1000, 1500))
abline(v = pick.peaks(lcms[1,,1], 20), col = "blue")
} else {
  cat("Package ptw not available.\nInstall it by typing 'install.packages(\"ptw\")'")
}
}
\keyword{manip}
