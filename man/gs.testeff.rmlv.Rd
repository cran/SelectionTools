\name{gs.testeff.rmlv}
\alias{gs.testeff.rmlv}
\title{Test estimated marker effects using the indicated method}

\description{
Test estimated marker effects using the indicated method.
}

\usage{
gs.testeff.rmlv(nperm = 1000, maxiter = 1000, precision = 0.001, hsq = 0.9,
                out.filename = "test.rmlv", auxfiles = TRUE,
                data.set ="default")
}

\arguments{
   \item{nperm}{Function argument used by this operation.}
   \item{maxiter}{Maximum number of iterations.}
   \item{precision}{Convergence tolerance.}
   \item{hsq}{Heritability parameter.}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the RMLV effect-test results when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
