\name{gs.esteff.rr}
\alias{gs.esteff.rr}
\title{Estimate marker effects using the RR method}

\description{
Estimate marker effects using the RR method.
}

\usage{
gs.esteff.rr(method = "BLUP", maxiter = 1000, precision = 0.0001, hsq = 0.80,
             alpha = 1, out.filename = "eff.rr", auxfiles = FALSE,
             data.set = "default")
}

\arguments{
   \item{method}{Method to use.}
   \item{maxiter}{Maximum number of iterations.}
   \item{precision}{Convergence tolerance.}
   \item{hsq}{Heritability parameter.}
   \item{alpha}{Method-specific tuning or significance parameter.}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the estimated marker effects from ridge regression when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
