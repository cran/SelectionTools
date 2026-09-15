\name{genome.parameter.set}
\alias{genome.parameter.set}
\title{Set genome parameters used by the simulation code}

\description{
Defines the chromosome structure used by the simulation backend.
}

\usage{
genome.parameter.set(no.chrom, no.hom=2, chrom.len)
}

\arguments{
\item{no.chrom}{Positive number of chromosomes.}
   \item{no.hom}{Number of homologues; the default is two.}
   \item{chrom.len}{Numeric vector containing one chromosome length for each chromosome.}
}
\value{
An invisible implementation-level \code{list} returned by the compiled state-setting or information routine. The function is primarily called for its side effect of resetting populations/map as needed and setting the genome parameters.
}

\details{
The genome definition is simulation-wide and determines the chromosome and homologue structure expected by population initialization, meiosis, crossing, single-seed descent, doubled-haploid production, and linkage-map operations.

\code{chrom.len} supplies one length per chromosome. Operations that model ordinary diploid marker data use two homologues. Changing genome parameters should therefore be coordinated with the loaded simulation linkage map and existing simulation populations.
}
