\name{PCA}
\alias{PCA}
\alias{scores}
\alias{loadings}
\alias{scores.PCA}
\alias{loadings.PCA}
\alias{variances}
\alias{project}
\alias{reconstruct}
\alias{summary.PCA}
\title{Principal Component Analysis}
\description{
  Functions for PCA: creating a PCA object, extracting variances, scores and
  loadings for individual PCs, projecting new data in the PC space,
  and reconstruction using a limited number of PCs.
}

\usage{
PCA(X, warn = TRUE)
\method{summary}{PCA}(object, varperc = 90, pc.select = c(1:5,10), \dots)
variances(object, npc = maxpc)
\method{scores}{PCA}(object, npc = maxpc, \dots)
\method{loadings}{PCA}(object, npc = maxpc, \dots)
reconstruct(object, npc = maxpc)
project(object, npc = maxpc, newdata, ldngs)
}

\arguments{
  \item{X}{a matrix, with each row representing an object.}
  \item{warn}{logical, whether or not to give a warning when the data
    are not mean-centered.}
  \item{object}{an object of class "PCA" (see below).}
  \item{varperc}{variance threshold in the \code{summary} function.}
  \item{\dots}{extra arguments, e.g., for printing the variance table
  (digits = ...).}
  \item{pc.select}{PCs to be included in the \code{summary} function.}
  \item{npc}{the number of PCs to be returned.}
  \item{newdata}{data (with the same number of variables as the original
    data) that are to be projected into the space of the first
    \code{npc} PCs.}
  \item{ldngs}{loadings to be used; by default the PCA loadings.}
}
\value{
  Function \code{PCA} returns an object of class "PCA" with components
  \item{scores}{object weights per PC.}
  \item{loadings}{variable weights per PC.}
  \item{var}{variance explained per PC.}
  \item{totalvar}{The total variance in the data set.}

  Function \code{summary.PCA} gives a short summary of the PCA model,
  stating how many PCs are needed to cover a certain percentage of the
  total variance, and for selected PCs gives the (cumulative) variance
  explained.

  Function \code{variances} returns the variances associated with each PC.

  Function \code{scores} returns the scores associated with each PC.

  Function \code{loadings} returns the loadings associated with each PC.

  Function \code{reconstruct} returns the reconstruction of the original
  data matrix, based on \code{npc} PCs.

  Function \code{project} projects the new data into the subspace
  spanned by the given loadings. If argument \code{ldngs} is given,
  arguments \code{pcamod} and \code{npc} are not needed.
}
\author{Ron Wehrens}
\references{
  R. Wehrens. "Chemometrics with R - Multivariate Data Analysis in the
  Natural Sciences and Life Sciences". 2nd edition, Springer, Heidelberg, 2019.
}
\seealso{\code{\link{plot.PCA}}}
\examples{
data(wines, package = "kohonen")
wines.PC <- PCA(scale(wines))
}
\keyword{multivariate}
