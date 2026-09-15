\name{gs.testeff.lsq}
\alias{gs.testeff.lsq}
\title{Test estimated marker effects using the indicated method}

\description{
Test estimated marker effects using the indicated method.
}

\usage{
gs.testeff.lsq(nperm = 0, out.filename = "test.lsq", auxfiles = TRUE,
               data.set ="default")
}

\arguments{
   \item{nperm}{Function argument used by this operation.}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the LSQ effect-test results when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
