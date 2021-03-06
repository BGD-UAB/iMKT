% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/completeMKT.R
\name{completimputedMKT}
\alias{completimputedMKT}
\title{Complete MK methodologies}
\usage{
completimputedMKT(daf, divergence, xlow, xhigh, seed)
}
\arguments{
\item{daf}{data frame containing DAF, Pi and P0 values}

\item{divergence}{data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes}

\item{xlow}{lower limit for asymptotic alpha fit}

\item{xhigh}{higher limit for asymptotic alpha fit}

\item{seed}{seed value (optional). No seed by default}
}
\value{
List with the diverse MKT results: standardMKT, FWW, DGRP, asymptoticMKT, iMKT
}
\description{
MKT calculation using all methodologies included in the package: standardMKT, FWW, DGRP, asymptoticMKT, iMKT.
}
\details{
Perform all MKT derived methodologies at once using the same input data and parameters.
}
\examples{
completimputedMKT(myDafData, myDivergenceData, xlow=0, xhigh=0.9)

}
\keyword{MKT}
