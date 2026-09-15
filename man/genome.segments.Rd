\name{genome.segments}
\alias{genome.segments}
\title{Extract or summarize genome segments for populations}

\description{
Extract or summarize genome segments for populations.
}

\usage{
genome.segments(pops, allele, chromosome=NULL, begin=NULL, end=NULL)
}

\arguments{
   \item{pops}{Population name or names.}
   \item{allele}{Allele specification.}
   \item{chromosome}{Chromosome identifier.}
   \item{begin}{Beginning position of a region.}
   \item{end}{Ending position of a region.}
}
\value{
A numeric \code{matrix} with one row per population and columns \code{Obs}, \code{Mean}, \code{SDev}, \code{Min}, \code{Q10}, \code{Med}, and \code{Max}, summarizing the number of genome segments.
}
