\name{genome.parameter.get}
\alias{genome.parameter.get}
\title{Returns the current genome parameters}

\description{
Returns the chromosome structure currently installed in the simulation backend.
}

\usage{
genome.parameter.get()
}

\details{
The result describes the simulation-wide genome parameters used to interpret chromosome segments and linkage-map positions. It includes the number of chromosomes, the number of homologues, and the chromosome-length vector.
}

\value{
A named \code{list} with components \code{no.chrom} (number of chromosomes), \code{no.hom} (number of homologues), and \code{chrom.len} (chromosome-length vector).
}
