\name{v-tropmaize-pop}
\alias{v-tropmaize-pop}
\alias{v.tropmaize.pop}
\docType{data}
\title{Genotype Data for Tropical Maize}
\description{
Compact, lossless representation of genotype data for 1135 markers and
264 tropical maize individuals.
}
\usage{
data("v-tropmaize-pop")
}
\format{
A list with 4 components:
\describe{
  \item{\code{marker}}{Character vector of 1135 marker names.}
  \item{\code{individual}}{Character vector of 264 individual identifiers.}
  \item{\code{genotype.levels}}{Character vector giving the genotype text
  represented by raw codes 0, 1, and 2, respectively: \code{"-1/-1"},
  \code{"1/1"}, and \code{"2/2"}.}
  \item{\code{genotype}}{A 1135 by 264 raw matrix. Each element is a raw
  code in \code{as.raw(0:2)} and indexes \code{genotype.levels} after
  conversion to integer and addition of one. Rows correspond to markers
  and columns correspond to individuals.}
}
}
\source{
Prepared from the original \file{ex-crossa.pop} text file.
}
\keyword{datasets}
