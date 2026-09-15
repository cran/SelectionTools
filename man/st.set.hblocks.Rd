\name{st.set.hblocks}
\alias{st.set.hblocks}
\title{Set haplotype-block information}

\description{
Set haplotype-block information.
}

\usage{
st.set.hblocks(haplotype.list, hap.symbol = "s", out.filename = "blocks",
               auxfiles = TRUE, data.set = "default")
}

\arguments{
   \item{haplotype.list}{Function argument used by this operation.}
   \item{hap.symbol}{Function argument used by this operation.}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{Pos}, \code{Name}, \code{Class}, and \code{Markers} describing the installed haplotype blocks when \code{auxfiles = TRUE} and the operation succeeds; otherwise \code{NULL}.
}
