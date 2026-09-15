\name{evaluate.allele}
\alias{evaluate.allele}
\title{evaluate.allele()}

\description{
  Calculates allele frequencies
}

\usage{
evaluate.allele(PopName, Loci, missing=0)
}

\arguments{
   \item{PopName}{Name of the population to be evaluated}
   \item{Loci}{Names of the loci to be evaluated}
   \item{missing}{Function argument.}
}
\value{
A \code{data.frame} containing one column for each evaluated locus plus a \code{count} column giving the number of observations for each allele combination. Returns \code{NULL} invisibly if no evaluable loci/result are available.
}

\note{
 Only positive alleles (without zero) are allowed!
}
