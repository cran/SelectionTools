\name{gs.estimate.gv}
\alias{gs.estimate.gv}
\title{Estimate genomic values for an estimation or prediction data set}

\description{
Estimate genomic values for an estimation or prediction data set.
}

\usage{
gs.estimate.gv(estimation.set ="default", prediction.set ="default",
               out.filename = "breeding.values", auxfiles = TRUE)
}

\arguments{
   \item{estimation.set}{Name of the estimation data set.}
   \item{prediction.set}{Name of the prediction data set.}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
}
\value{
Invisibly returns a \code{data.frame} with genomic-value results when \code{auxfiles = TRUE} and estimation succeeds. It contains individual identifiers and predicted genomic values, and also observed phenotypes when those are present in the generated result. Otherwise \code{NULL}.
}
