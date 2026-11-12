\name{gs.build.V}
\alias{gs.build.V}
\title{Builds the genomic-selection variance or covariance matrix for a data set}

\description{
Builds the genomic-selection variance or covariance matrix for a data set.
}

\usage{
gs.build.V(out.filename = "V.matrix", auxfiles = FALSE, data.set = "default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the generated variance/covariance matrix \code{V} when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
