\name{shootout}
\alias{shootout}
\docType{data}
\title{Shootout NIR data}
\description{NIR data from 654 tablets, measured at two different
  instruments. The data have been divided in training, test and
  validation sets.
}
\usage{data(shootout)}
\format{Variable \code{shootout} is a list containing spectral data of
  tablets, measured on two instruments, as well as response variables.
}
\details{For every tablet, three reponse variables are measured: the
  amount of active ingredient (nominally 200 mg/tablet), the weight and
  the hardness. Typically, one wants to estimate the amount of active
  ingredient from the NIR spectra, a straightforward multivariate
  calibration problem. The goal of the shootout competition was to find
  the optimal way to transfer a calibration model of the first
  instrument to the second.
}
\references{
  R. Wehrens. "Chemometrics with R - Multivariate Data Analysis in the
  Natural Sciences and Life Sciences". 2nd edition, Springer, Heidelberg, 2019.
}
\source{http://www.idrc-chambersburg.org/shootout2002.html}
\examples{
data(shootout)
plot(seq(600, 1898, by = 2), shootout$calibrate.1[1,], type = "l",
xlab = "wavelength", ylab = "log(1/T)")
}
\keyword{datasets}
