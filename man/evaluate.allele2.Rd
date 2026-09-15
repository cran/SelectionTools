\name{evaluate.allele2}
\alias{evaluate.allele2}
\title{Evaluate allele information for selected loci in a population}

\description{
Evaluate allele information for selected loci in a population.
}

\usage{
evaluate.allele2(PopName, Loci, missing=0)
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{Loci}{Locus or loci to process.}
   \item{missing}{Value or code used for missing data.}
}
\value{
A \code{data.frame} containing one column for each evaluated locus plus a \code{count} column giving the number of observations for each allele combination. Returns \code{NULL} invisibly if no evaluable loci/result are available.
}
