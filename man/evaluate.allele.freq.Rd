\name{evaluate.allele.freq}
\alias{evaluate.allele.freq}
\title{Calculate or return allele frequencies for selected loci}

\description{
Calculate or return allele frequencies for selected loci.
}

\usage{
evaluate.allele.freq(PopName, Loci, missing=0)
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{Loci}{Locus or loci to process.}
   \item{missing}{Value or code used for missing data.}
}
\value{
A \code{data.frame} containing the allele counts returned by \code{evaluate.allele()} and an added \code{freq} column containing relative frequencies. Returns \code{NULL} invisibly when the request cannot be evaluated.
}
