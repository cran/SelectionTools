\name{gd.pcoa}

\alias{gd.pcoa}

\title{gd.pcoa()}

\description{Calculates a principal coordinate analysis for marker
  data. It is a wrapper for the \code{cmdscale} function. }

\usage{
gd.pcoa(dta, k=3)
}

\arguments{
   \item{dta}{Name of the dataset containing the observed marker
   genotypes}
   \item{k}{Number of the principal coordinates that should be returned}
}

\value{
A numeric \code{matrix} of principal-coordinate scores returned by \code{cmdscale()}, with rows corresponding to populations and columns to the first \code{k} principal coordinates.
}
