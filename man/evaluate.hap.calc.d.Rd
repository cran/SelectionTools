\name{evaluate.hap.calc.d}
\alias{evaluate.hap.calc.d}
\title{Calculate a haplotype statistic for a pair of loci}

\description{
Calculate a haplotype statistic for a pair of loci.
}

\usage{
evaluate.hap.calc.d(locusA, locusB)
}

\arguments{
   \item{locusA}{First locus.}
   \item{locusB}{Second locus.}
}
\value{
An invisible numeric vector of length three, in the order \code{dd}, \code{dp}, and \code{rr}, containing the values returned by the compiled haplotype-statistic calculation.
}
