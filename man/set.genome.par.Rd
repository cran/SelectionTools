\name{set.genome.par}
\alias{set.genome.par}
\title{set.genome.par()}

\description{
Defines simulation chromosome number, ploidy, and chromosome lengths.
}

\usage{
set.genome.par(no.chrom, no.hom, chrom.len)
}

\arguments{
\item{no.chrom}{Positive number of simulation chromosomes.}
   \item{no.hom}{Number of homologues. Diploid population operations use two.}
   \item{chrom.len}{Numeric vector of chromosome lengths with one value per chromosome.}
}
\value{
Returns the same implementation-level \code{list} as \code{genome.parameter.set()}; the main effect is setting the simulation genome parameters.
}

\details{
The genome definition controls the chromosome structure used by simulation-population operations. \code{no.chrom} gives the number of chromosomes, \code{no.hom} the number of homologues, and \code{chrom.len} the chromosome lengths in the simulation map unit.

Population-generation operations including crossing, single-seed descent, doubled haploids, and the SelectionTools marker-to-simulation bridge require a diploid genome with two homologues. The chromosome-length vector must correspond to the declared chromosome count.

This function and \code{genome.parameter.set} modify the same simulation genome state.
}
