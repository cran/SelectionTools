\name{gd.data.parameters}

\alias{gd.data.parameters}

\title{gd.data.parameters()}

\description{
  Determines the structure of a data set containing molecular marker data.
}

\usage{
gd.data.parameters(dta)
}

\arguments{
   \item{dta}{Name of the dataset contaning the observed marker
   genotypes. The format of the dataset is illustrated by the datasets
   gd.dummy.populations and gd.maize.populations.}
}

\value{
A named \code{list} describing the genetic-diversity data object, with components \code{colhead}, \code{rowhead}, \code{no.mar}, \code{no.all}, \code{no.pop}, \code{no.ind}, \code{ind.list}, \code{mar.list}, and \code{pop.list}. These contain parsed column and row headings, marker/allele/population counts, individual indices, marker names, and population names, respectively.
}
