\name{bdata}
\alias{bdata}
\docType{data}
\title{HPLC-UV data of two chemical mixtures}
\description{Two chemical mixtures of three compounds have been measured
  using HPLC-UV. Two of the compounds are known: diazinon and
  parthion-ethyl, both organophosphorus pesticides. Each data matrix
  consists of 73 wavelengths and 40 time points. The challenge is to
  infer the pure spectra of the individual compounds, as well as their
  time profiles.
}
\usage{data(bdata)}
\format{A list of four elements. The first two, \code{d1} and \code{d2},
  are the mixture matrices of the two analytes and one unknown
  interferent. The last two, \code{sp1} and \code{sp2}, contain the pure
  spectra of the two analytes.
}
\source{
Original matlab data files obtained from
\code{http://www.ub.edu/mcr/web_mcr/download_dataHPLC.html} (bdataset.zip).
No longer available.
}
\references{
  R. Wehrens. "Chemometrics with R - Multivariate Data Analysis in the
  Natural Sciences and Life Sciences". Springer, 2nd edition, Heidelberg, 2019.

  R. Tauler, S. Lacorte and D. Barcelo. "Application of multivariate curve
self-modeling curve resolution for the quantitation of trace levels of
organophosphorous pesticides in natural waters from interlaboratory
studies". J. of Chromatogr. A, 730, 177-183 (1996).
}
\examples{
data(bdata)
persp(bdata$d1, phi = 20, theta = 34, expand = .5,
      xlab = "Time", ylab = "Wavelength")
}
\keyword{datasets}
