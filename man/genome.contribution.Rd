\name{genome.contribution}
\alias{genome.contribution}
\title{genome.contribution()}

\description{
 Calculates the distribution of one allele in populations.
}

\usage{
genome.contribution(pops, allele, chromosome=NULL, begin=NULL, end=NULL)
}

\arguments{
   \item{pops}{Names of the populations separated by blanks}
   \item{allele}{Allele which is considered}
   \item{chromosome}{Chromosome}
   \item{begin}{Start of chromosome region to investigate}
   \item{end}{End of the chromosome region}
}
\value{
A numeric \code{matrix} with one row per population and columns \code{Obs}, \code{Mean}, \code{SDev}, \code{Min}, \code{Q10}, \code{Med}, and \code{Max}, summarizing genome-contribution values.
}

\details{
Analysis is done on the level of the chromosomes, independently of the
loaded marker map. Summarizes the results of \code{evaluate.genome}.
}

\note{Further details appear in the \code{lib00.example1} for an application of the command.}
