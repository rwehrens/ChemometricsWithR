\name{Prostate2000Raw}
\alias{Prostate2000Raw}
\title{Prostate Cancer 2000 Raw Spectra}
\concept{prostate cancer}
\description{A data object of class \code{msSet},
consisting of 654 mass spectra (327 spectra in duplicate) from 2000 to
20000 Da, which were generated from patients with prostate cancer,
benign prostatic hypertrophy, and normal controls. These spectra are
already baseline corrected and normalized. Please see the references for
more details.

Since the original package msProstate is orphaned at the end of 2012,
the data are included in the ChemometricsWithR package so that the
examples in the book are still executable. This manual page has been
adapted to reflect this.
} 

\references{
B.L. Adam, Y. Qu, J.W. Davis, M.D. Ward, M.A. Clements, L.H. Cazares,
O.J. Semmes, P.F. Schellhammer, Y. Yasui, Z. Feng, and G.L. Wright, Jr.,
"Serum protein fingerprinting coupled with a pattern-matching algorithm
distinguishes prostate cancer from benign prostate hyperplasia and
healthy men," \emph{Cancer Research}, 62(13):3609--14, 2002.

Y. Qu, B.L. Adam, Y. Yasui, M.D. Ward, L.H. Cazares,
P.F. Schellhammer, Z. Feng, O.J. Semmes, and G.L. Wright Jr.,
"Boosted decision tree analysis of surface-enhanced laser
desorption/ionization mass spectral serum profiles discriminates
prostate cancer from noncancer patients",
\emph{Clinical Chemistry}, 48(10):1835--43, 2002.

 R. Wehrens, "Chemometrics with R - Multivariate Data Analysis in the
  Natural Sciences and Life Sciences". 2nd edition, Springer, Heidelberg, 2019.
}

\examples{
## Examples have been changed from the original man page upon inclusion
## in the ChemometricsWithRData package
data("Prostate2000Raw")

## plot a few spectra, partially
matplot(Prostate2000Raw$mz[1:8000],
        Prostate2000Raw$intensity[1:8000, 1:5], type = "l",
        lty = 1, col = 1:5, xlab = "m/z", ylab = "response")
}
\keyword{datasets}
