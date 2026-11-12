\name{st.set.hblocks}
\alias{st.set.hblocks}
\title{Sets haplotype-block information}

\description{
Sets haplotype-block information.
}

\usage{
st.set.hblocks(haplotype.list, hap.symbol = "s", out.filename = "blocks",
               auxfiles = TRUE, data.set = "default")
}

\arguments{
   \item{haplotype.list}{Function argument used by this operation.}
   \item{hap.symbol}{Function argument used by this operation.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{Pos}, \code{Name}, \code{Class}, and \code{Markers} describing the installed haplotype blocks when \code{auxfiles = TRUE} and the operation succeeds; otherwise \code{NULL}.
}
