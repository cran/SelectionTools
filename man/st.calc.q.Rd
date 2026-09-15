\name{st.calc.q}
\alias{st.calc.q}
\title{Calculate the package-specific Q statistic or quantity}

\description{
Calculate the package-specific Q statistic or quantity.
}

\usage{
st.calc.q(pop.type = "DH", t = 0, auxfiles = TRUE, data.set ="default")
}

\arguments{
   \item{pop.type}{Population type.}
   \item{t}{Generation or method-specific time parameter.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{Locus1}, \code{Locus2}, and \code{q}, containing the calculated pairwise \code{q} values; otherwise \code{NULL}.
}
