\name{evaluate.all.loci}
\alias{evaluate.all.loci}
\title{Evaluates all loci in a population}

\description{
Evaluates all loci in a population.
}

\usage{
evaluate.all.loci(PopName,
                  loci = paste(as.vector(unlist(linkage.map.get()[,3])), collapse=" "))
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{loci}{Locus or loci to process.}
}
\value{
The same \code{data.frame} structure as \code{evaluate.locus()}, evaluated for all loci in the current map: rows are loci and columns are \code{homozygote}, \code{heterozygote}, and \code{missing}.
}
