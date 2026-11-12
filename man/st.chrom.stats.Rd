\name{st.chrom.stats}
\alias{st.chrom.stats}
\title{Returns chromosome statistics for a data set}

\description{
Returns chromosome statistics for a data set.
}

\usage{
st.chrom.stats(filename = "chrom.stats", auxfiles = TRUE, data.set ="default")
}

\arguments{
   \item{filename}{Character string retained for compatibility. Results are returned from memory without creating a file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{NLoci}, and \code{Length}, giving chromosome identifier, number of loci, and chromosome length; otherwise \code{NULL} if the calculation fails.
}
