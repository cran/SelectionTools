\name{evaluate.locus}
\alias{evaluate.locus}
\title{Evaluate selected loci in a population}

\description{
Evaluate selected loci in a population.
}

\usage{
evaluate.locus(PopName, eLoci)
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{eLoci}{Locus or loci to evaluate.}
}
\value{
A \code{data.frame} with one row for each requested locus and columns \code{homozygote}, \code{heterozygote}, and \code{missing}, containing the corresponding genotype counts.
}
