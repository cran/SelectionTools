\name{evaluate.genotype2}
\alias{evaluate.genotype2}
\title{Evaluate genotype information for selected loci in a population}

\description{
Evaluate genotype information for selected loci in a population.
}

\usage{
evaluate.genotype2(PopName, Loci, alleles=0, mode=2, bg=0)
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{Loci}{Locus or loci to process.}
   \item{alleles}{Allele specification or coding.}
   \item{mode}{Evaluation or operation mode.}
   \item{bg}{Background or missing-data code used by the evaluator.}
}
\value{
A \code{data.frame} containing genotype combinations and the count/frequency information returned by the compiled genotype evaluator. Returns \code{NULL} invisibly if no result is available.
}
