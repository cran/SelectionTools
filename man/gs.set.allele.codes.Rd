\name{gs.set.allele.codes}
\alias{gs.set.allele.codes}
\title{Set allele codes used by the genomic-selection routines}

\description{
Set allele codes used by the genomic-selection routines.
}

\usage{
gs.set.allele.codes(aa = 0, aA = 1, AA = 2, an = 0.5, nn = 1, nA = 1.5,
                    data.set ="default")
}

\arguments{
   \item{aa}{Function argument used by this operation.}
   \item{aA}{Function argument used by this operation.}
   \item{AA}{Function argument used by this operation.}
   \item{an}{Function argument used by this operation.}
   \item{nn}{Function argument used by this operation.}
   \item{nA}{Function argument used by this operation.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
\code{NULL}, returned invisibly. The function is called for its side effect of setting the allele codes used by the genomic-selection routines.
}
