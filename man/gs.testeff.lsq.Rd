\name{gs.testeff.lsq}
\alias{gs.testeff.lsq}
\title{Tests estimated marker effects using the indicated method}

\description{
Tests estimated marker effects using the indicated method.
}

\usage{
gs.testeff.lsq(nperm = 0, out.filename = "test.lsq", auxfiles = TRUE,
               data.set ="default")
}

\arguments{
   \item{nperm}{Function argument used by this operation.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the LSQ effect-test results when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
