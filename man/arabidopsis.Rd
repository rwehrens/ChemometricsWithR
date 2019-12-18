\name{arabidopsis}
\alias{arabidopsis}
\alias{arabidopsis.Y}
\docType{data}
\title{LC-MS metabolomics data sets from Arabidopsis samples}
\description{LC-MS data from Arabidopsis thaliana samples. The
  \code{arabidopsis} data object contains relative 
  intensities of 567 reconstructed metabolites (columns) for 761
  samples (rows). A sizeable fraction of intensities are missing, in
  most cases because the corresponding metabolites are below the
  detection level. The corresponding meta-information object
  (\code{arabidopsis.Y})contains for every sample batch and 
  sequence information, as well as (coded) information on the genotype
  and the sample type (study sample or reference sample). Processing of
  the raw data has been done with Metalign and MSClust programs.
}
\usage{data(arabidopsis)}
\references{
  "@Article{Wehrens2016,
  author =	 {Ron Wehrens and Jos.~A.~Hageman and Fred~van~Eeuwijk
                  and Rik~Kooke and P\'adraic~J.~Flood and Erik
                  Wijnker and Joost~J.B.~Keurentjes and Arjen~Lommen
                  and Henri\"ette~D.L.M.~van~Eekelen and
                  Robert~D.~Hall and Roland~Mumm and Ric~C.H.~de~Vos},
  title =	 {Improved batch correction in untargeted MS-based
                  metabolomics},
  journal =	 {Metabolomics},
  year =	 2016,
  volume =	 12,
  DOI =		 {10.1007/s11306-016-1015-8},
  pages =	 {1--12}
}"
}
\examples{
data(arabidopsis)
dim(arabidopsis)
sum(is.na(arabidopsis))
dim(arabidopsis.Y)
}
\keyword{datasets}
