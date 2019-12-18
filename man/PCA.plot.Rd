\name{PCA.plot}
\alias{plot.PCA}
\alias{scoreplot}
\alias{loadingplot}
\alias{scoreplot.PCA}
\alias{loadingplot.PCA}
\alias{screeplot}
\alias{biplot.PCA}
\title{Principal Component Analysis plotting functions}
\description{
  Plotting functions for PCA: for scores, loadings, scores and loadings
  simultaneously (a biplot), and variances (a screeplot, where the log
  of the explained variance is plotted for each PC).
}

\usage{
\method{scoreplot}{PCA}(object, pc = c(1, 2), pcscores = scores(object),
          show.names = FALSE, xlab, ylab, xlim, ylim, \dots)
\method{loadingplot}{PCA}(object, pc = c(1, 2), pcloadings = loadings(object),
            scalefactor = 1, add = FALSE, show.names = FALSE,
            xlab, ylab, xlim, ylim, col = "blue", min.length =
            0.01, varnames = NULL, \dots)
\method{biplot}{PCA}(x, pc = c(1,2),
        show.names = c("none", "scores", "loadings", "both"),
        score.col = 1, loading.col = "blue",
        min.length = .01, varnames = NULL, \dots)
screeplot(object, type = c("scree", "percentage"), npc, \dots)
}

\arguments{
  \item{x, object}{an object of class "PCA" (see below).}
  \item{pc}{which PCs to show.}
  \item{pcscores}{matrix of scores, by default the scores of the PCA
  model object.}
  \item{show.names}{show names rather than plotting
    symbols. For \code{loadingplot} and \code{scoreplot} a logical
    (default: FALSE), for \code{biplot} one of 'scores', 'loadings',
    'both' or 'none' (default).}
  \item{xlab, ylab, xlim, ylim, col}{graphical parameters of the plot.}
  \item{pcloadings}{matrix of loadings, by default the loadings of the
    PCA model object.}
  \item{scalefactor}{scaling factor for the loadings; used internally, when
    the \code{loadingplot} function is called from within
    \code{biplot.PCA}.}
  \item{add}{logical, whether to add to the existing plot (again, useful
    when \code{loadingplot} is called from within \code{biplot.PCA}).}
  \item{npc}{how many PCs to show in the scree plot (starting from 1).}
  \item{type}{show a real screeplot (\code{scree}) or show the
    percentage of variance explained (\code{percentage}).}
  \item{score.col, loading.col}{colours of the scores and loadings in a
    biplot.}
  \item{min.length}{minimal length of loading vectors to be plotted by
    arrows. Vectors that are too short lead to warning messages, are not
    interesting, and only clutter the graphic.}
  \item{varnames}{alternative vector of variable names.}
  \item{\dots}{Graphical arguments passed on to lower-level plotting
  functions.}
}
\details{Score plots and loading plots show the amount of explained
  variance at the axis labels only when PCA
  has been performed at mean-centered data.}
\author{Ron Wehrens}
\references{
  R. Wehrens. "Chemometrics with R - Multivariate Data Analysis in the
  Natural Sciences and Life Sciences". 2nd edition, Springer, Heidelberg, 2019.
}
\seealso{\code{\link{PCA}}}
\examples{
data(wines, package = "kohonen")
wines.PC <- PCA(scale(wines))
wine.classes <- as.integer(vintages)
scoreplot(wines.PC, col = wine.classes, pch = wine.classes)
loadingplot(wines.PC, show.names = TRUE)
biplot(wines.PC, score.col = wine.classes)
screeplot(wines.PC)
}
\keyword{hplot}
