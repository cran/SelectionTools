\name{gs.estimate.gv}
\alias{gs.estimate.gv}
\title{Estimates genomic values for an estimation or prediction data set}

\description{
Estimates genomic values for an estimation or prediction data set.
}

\usage{
gs.estimate.gv(estimation.set ="default", prediction.set ="default",
               out.filename = "breeding.values", auxfiles = TRUE)
}

\arguments{
   \item{estimation.set}{Name of the estimation data set.}
   \item{prediction.set}{Name of the prediction data set.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
}
\value{
Invisibly returns a \code{data.frame} with genomic-value results when \code{auxfiles = TRUE} and estimation succeeds. It contains individual identifiers and predicted genomic values, and also observed phenotypes when those are present in the generated result. Otherwise \code{NULL}.
}
