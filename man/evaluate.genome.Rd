\name{evaluate.genome}
\alias{evaluate.genome}
\title{evaluate.genome()}

\description{
  Computes the frequency
  distribution of an allele across the whole genome.
}

\usage{
evaluate.genome(PopName, allele, chromosome=NULL, begin=NULL, end=NULL,
                verbose=FALSE)
}

\arguments{
   \item{PopName}{Name of the population to be evaluated}
   \item{allele}{Function argument.}
   \item{chromosome}{Chromosome}
   \item{begin}{Start of chromosome region to investigate}
   \item{end}{End of the chromosome region}
   \item{verbose}{Function argument.}
}

\details{
Analysis is done on the level of the chromosomes, independently of the
loaded marker map. 
}


\value{
A \code{data.frame} with one row per individual. It always contains \code{ratio}, the evaluated allele proportion. With \code{verbose = TRUE} it additionally contains \code{length}, \code{total}, and \code{noBlocks}, as returned by the genome-evaluation routine.
}

\note{Further details appear in the \code{lib00.example1} for an application of the command.}
