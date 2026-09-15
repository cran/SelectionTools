\name{st.calc.ld}
\alias{st.calc.ld}
\title{Calculate linkage disequilibrium from marker data}

\description{
Calculate linkage disequilibrium from marker data.
}

\usage{
st.calc.ld(ld.measure = "r2", auxfiles = TRUE, data.set = "default")
}

\arguments{
   \item{ld.measure}{Function argument used by this operation.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{Locus1}, \code{Locus2}, \code{Name1}, \code{Name2}, and \code{LD}, containing pairwise linkage-disequilibrium results; otherwise \code{NULL}.
}
