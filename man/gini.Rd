\name{gini}
\alias{gini}
\title{Gini impurity index for cart objects}
\description{Calculation of the change in the Gini impurity index for
  classification and regression trees. The function returns changes in
  the gini index associated with using individual values of \code{x} as
  split points. Included for demonstration purposes.
}
\usage{
gini(x, clss)
}
\arguments{
  \item{x}{Numeric vector of length n.}
  \item{clss}{Clss labels, length n.}
}
\references{
  R. Wehrens. "Chemometrics with R - Multivariate Data Analysis in the
  Natural Sciences and Life Sciences". 2nd edition, Springer, Heidelberg, 2019.
}
\value{The change in Gini impurity index, given a vector of 
  possible splits, and a vector of clss labels. Lower values indicate
  more pure leaves.}
\author{Ron Wehrens}
\keyword{manip}
