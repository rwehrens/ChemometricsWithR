\name{MCR}
\alias{MCR}
\alias{mcr}
\alias{efa}
\title{Functions for Multivariate Curve Resolution}
\usage{
mcr(x, init, what = c("row", "col"), convergence = 1e-08,
    maxit = 50)
efa(x, ncomp)
}
\description{Multivariate Curve Resolution, or MCR, decomposes a
  bilinear matrix into its pure components. A classical example is a
  matrix consisting of a series of spectral measurements on a mixture of
  chemicals for following the reaction. At every time point, a spectrum
  is measured that is a linear combination of the pure spectra. The goal
  of MCR is to resolve the pure spectra and concentration profiles over
  time.}
\arguments{
  \item{x}{Data matrix}
  \item{init}{Initial guess for pure compounds}
  \item{what}{Whether the pure compounds are rows or columns of the data
    matrix}
  \item{convergence}{Convergence criterion}
  \item{maxit}{Maximal number of iterations}
  \item{ncomp}{Number of pure compounds}
}
\details{
MCR uses repeated application of least-squares regression to find pure
profiles and spectra. The method is iterative; EFA is a
method to provide initial guesses.
}
\value{Function \code{mcr} returns a list containing
  \item{C}{An estimate of the pure "concentration profiles"}
  \item{S}{An estimate of the pure "spectra"}
  \item{resids}{The residuals of the final decomposition}
  \item{rms}{Root-mean-square values of the individual iterations}
  
  Function \code{efa} returns a list containing
  \item{pure.compounds: }{A matrix containing \code{ncomp} pure
    compounds, usually concentration profiles at specific wavelengths}
  \item{forward: }{The development of the singular values of the
    reduced data matrix when increasing the number of columns in the
    forward direction}
  \item{backward: }{The development of the singular values of the
    reduced data matrix when increasing the number of columns in the
    backwarddirection}
}
\author{Ron Wehrens}
\references{
  R. Wehrens. "Chemometrics with R - Multivariate Data Analysis in the
  Natural Sciences and Life Sciences". 2nd edition, Springer, Heidelberg, 2019.
}
\examples{
data(bdata)
D1.efa <- efa(bdata$d1, 3)
matplot(D1.efa$forward, type = "l")
matplot(D1.efa$backward, type = "l")
matplot(D1.efa$pure.comp, type = "l")

D1.mcr.efa <- mcr(bdata$d1, D1.efa$pure.comp, what = "col")
matplot(D1.mcr.efa$C, type = "l", main = "Concentration profiles")
matplot(t(D1.mcr.efa$S), type = "l", main = "Pure spectra")
}
\keyword{manip}
