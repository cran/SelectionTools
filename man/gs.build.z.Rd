\name{gs.build.Z}
\alias{gs.build.Z}
\title{gs.build.Z()}
\description{
  The function \emph{gs.build.Z} builds the design (Z) matrix out of the marker data
  for genomic selection.
}
\usage{
gs.build.Z(out.filename = "Z.matrix", auxfiles = FALSE, data.set = "default")
}
\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{character. Name of the data set that is used by the function}
}

\details{
By default, genotypes are coded as \code{0}, \code{1}, or \code{2},
representing the number of copies of the effect allele. When effect
alleles have not already been defined, the larger observed allele
code at each marker is used as the effect allele. For biallelic
STvcf data with both REF and ALT observed, this corresponds to
\code{0} for REF/REF, \code{1} for REF/ALT, and \code{2} for ALT/ALT. 
Missing genotypes receive the same value as heterozygotes
(\code{1} by default), so these cannot be distinguished in the
resulting matrix.
}

\value{
Invisibly returns a \code{matrix} containing the generated marker/design matrix \code{Z} when \code{auxfiles = TRUE} and generation succeeds; otherwise \code{NULL}.
}


