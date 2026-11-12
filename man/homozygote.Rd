\name{homozygote}
\alias{homozygote}
\title{homozygote()}

\description{
  Produces a dataframe of a population
  whith \emph{NoInd} individuals which carry at all loci the same allele
  (\emph{allele}). This dataframe can be used to initialise a population.
}

\usage{
homozygote(allele, NoInd = 1)
}

\arguments{
   \item{allele}{Allele which is carried at all loci}
   \item{NoInd}{Number of individuals to be generated}
}
\value{
A \code{data.frame} with columns \code{ind}, \code{chrom}, \code{hom}, \code{pos}, and \code{all}, representing the requested homozygous genome.
}

\note{
  Further details appear in \code{init.population} for an example
}
