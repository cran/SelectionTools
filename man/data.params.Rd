\name{data.params}
\alias{data.params}
\title{Returns or inspect parameters associated with a data object}

\description{
Returns or inspect parameters associated with a data object.
}

\usage{
data.params(dta)
}

\arguments{
   \item{dta}{Input data object.}
}
\value{
A named \code{list} describing the supplied data object, with components:

  \item{colhead}{Parsed column headings identifying populations and individuals.}
  \item{rowhead}{Parsed row headings identifying markers and alleles.}
  \item{no.mar}{Number of markers.}
  \item{no.all}{Number of alleles per marker.}
  \item{no.pop}{Number of populations.}
  \item{no.ind}{Numbers of individuals per population.}
  \item{ind.list}{Individual numbers within populations.}
  \item{mar.list}{Marker names.}
  \item{pop.list}{Population names.}

}
