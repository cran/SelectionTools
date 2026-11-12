\name{evaluate.ld}
\alias{evaluate.ld}
\title{Evaluates linkage disequilibrium for a population and genomic position}

\description{
Evaluates linkage disequilibrium for a population and genomic position.
}

\usage{
evaluate.ld(PopName, allele, chrom, pos)
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{allele}{Allele specification.}
   \item{chrom}{Chromosome identifier.}
   \item{pos}{Position on the chromosome.}
}
\value{
A one-column \code{data.frame} named \code{length}, containing the linkage-disequilibrium segment lengths returned for the requested population, allele and position.
}
