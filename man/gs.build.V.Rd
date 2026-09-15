\name{gs.build.V}
\alias{gs.build.V}
\title{Build the genomic-selection variance or covariance matrix for a data set}

\description{
Build the genomic-selection variance or covariance matrix for a data set.
}

\usage{
gs.build.V(out.filename = "V.matrix", auxfiles = FALSE, data.set = "default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the generated variance/covariance matrix \code{V} when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
