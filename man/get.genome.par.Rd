\name{get.genome.par}
\alias{get.genome.par}
\title{get.genome.par()}

\description{
Returns the current simulation genome parameters.
}


\usage{
get.genome.par()
}

\details{
The returned values describe the genome used by simulation populations: chromosome count, number of homologues, and chromosome lengths. These parameters are separate from a SelectionTools marker data set, although the marker-to-simulation bridge can derive and install compatible simulation genome parameters from a marker map.

}

\value{
Returns the same named \code{list} as \code{genome.parameter.get()}, with components \code{no.chrom}, \code{no.hom}, and \code{chrom.len}.
}
