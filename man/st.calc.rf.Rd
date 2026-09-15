\name{st.calc.rf}
\alias{st.calc.rf}
\title{Calculate recombination fractions}

\description{
Calculate recombination fractions.
}

\usage{
st.calc.rf(map.function = "Haldane", auxfiles = TRUE, data.set = "default")
}

\arguments{
   \item{map.function}{Map function to use.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{Locus1}, \code{Locus2}, and \code{RecFreq}, containing pairwise recombination fractions; otherwise \code{NULL}.
}
