\name{evaluate.hap.calc}
\alias{evaluate.hap.calc}
\title{Calculate haplotype information for a pair of loci}

\description{
Calculate haplotype information for a pair of loci.
}

\usage{
evaluate.hap.calc(locusA, locusB)
}

\arguments{
   \item{locusA}{First locus.}
   \item{locusB}{Second locus.}
}
\value{
A numeric vector of length five, in the order \code{noTI}, \code{noI}, \code{nodH}, \code{dLocA}, and \code{dLocB}, containing the values returned by the compiled pairwise haplotype calculation.
}
