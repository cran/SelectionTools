\name{gs.return.pvals}
\alias{gs.return.pvals}
\title{Returns marker p-values for a genomic-selection data set}

\description{
Returns marker p-values for a genomic-selection data set.
}

\usage{
gs.return.pvals(out.filename ="effects", data.set ="default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
A \code{data.frame} containing the marker p-value table. If the underlying routine reports an error, \code{NULL} is returned invisibly. The intermediate file used to transfer the table is deleted after reading.
}
