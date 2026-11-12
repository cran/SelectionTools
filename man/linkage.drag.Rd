\name{linkage.drag}
\alias{linkage.drag}
\title{Calculates or summarize linkage drag}

\description{
Calculates or summarize linkage drag.
}

\usage{
linkage.drag(pops, allele, chrom, pos)
}

\arguments{
   \item{pops}{Population name or names.}
   \item{allele}{Allele specification.}
   \item{chrom}{Chromosome identifier.}
   \item{pos}{Position on the chromosome.}
}
\value{
A numeric \code{matrix} with one row per population and columns \code{Obs}, \code{Mean}, \code{SDev}, \code{Min}, \code{Q10}, \code{Med}, and \code{Max}, summarizing linkage-drag segment lengths.
}
