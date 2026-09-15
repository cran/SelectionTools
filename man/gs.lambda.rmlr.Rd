\name{gs.lambda.rmlr}
\alias{gs.lambda.rmlr}
\title{Estimate or set genomic-selection lambda values using the rmlr method}

\description{
Estimate or set genomic-selection lambda values using the rmlr method.
}

\usage{
gs.lambda.rmlr(maxiter = 1000, precision = 0.001, hsq = 0.9,
               out.filename = "lambda", auxfiles = TRUE, data.set = "default")
}

\arguments{
   \item{maxiter}{Maximum number of iterations.}
   \item{precision}{Convergence tolerance.}
   \item{hsq}{Heritability parameter.}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the lambda results from the RMLR method when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
