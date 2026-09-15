\name{evaluate.genotype}
\alias{evaluate.genotype}
\title{evaluate.genotype()}

\description{
  Calculates genotype frequencies
}

\usage{
evaluate.genotype(PopName, Loci, alleles=0, mode=2, bg=0)
}

\arguments{
   \item{PopName}{Name of the population to be evaluated}
   \item{Loci}{A List of Loci}
   \item{alleles}{A list of different Alleles}
   \item{mode}{Considering gametic phase or not}
   \item{bg}{Function argument.}
}
\value{
A \code{data.frame} containing the genotype combinations and the count/frequency information returned by the compiled genotype evaluator. If a specific allele combination is requested the result can be reduced to the matching row. Returns \code{NULL} invisibly if no result is available.
}
