\name{st.chrom.stats}
\alias{st.chrom.stats}
\title{Return chromosome statistics for a data set}

\description{
Return chromosome statistics for a data set.
}

\usage{
st.chrom.stats(filename = "chrom.stats", auxfiles = TRUE, data.set ="default")
}

\arguments{
   \item{filename}{Character string retained for compatibility and used as part of the name of the internal session-temporary result file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{NLoci}, and \code{Length}, giving chromosome identifier, number of loci, and chromosome length; otherwise \code{NULL} if the calculation fails.
}
