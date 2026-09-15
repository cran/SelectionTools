\name{evaluate.hap.init}
\alias{evaluate.hap.init}
\title{Initialize haplotype evaluation for a population}

\description{
Initialize haplotype evaluation for a population.
}

\usage{
evaluate.hap.init(PopName, doubleHetero=0, missing=0)
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{doubleHetero}{Control for handling double heterozygotes.}
   \item{missing}{Value or code used for missing data.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of initializing haplotype evaluation state.
}
