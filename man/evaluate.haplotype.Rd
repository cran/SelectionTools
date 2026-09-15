\name{evaluate.haplotype}
\alias{evaluate.haplotype}
\title{Evaluate haplotypes for selected loci in a population}

\description{
Evaluate haplotypes for selected loci in a population.
}

\usage{
evaluate.haplotype(PopName, Loci, missing=0,
                   bg=plabsim$allele.missing.indicator)
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{Loci}{Locus or loci to process.}
   \item{missing}{Value or code used for missing data.}
   \item{bg}{Background or missing-data code used by the evaluator.}
}
\value{
A \code{data.frame} containing the haplotype-evaluation table returned for the requested population and loci.
}
