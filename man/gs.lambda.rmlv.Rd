\name{gs.lambda.rmlv}
\alias{gs.lambda.rmlv}
\title{Estimates or sets genomic-selection lambda values using the rmlv method}

\description{
Estimates or sets genomic-selection lambda values using the rmlv method.
}

\usage{
gs.lambda.rmlv(maxiter = 1000, precision = 0.001, out.filename = "lambda",
               auxfiles = TRUE, data.set = "default")
}

\arguments{
   \item{maxiter}{Maximum number of iterations.}
   \item{precision}{Convergence tolerance.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the lambda results from the RMLV method when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
