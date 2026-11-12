\name{gs.lambda.rmla.02}
\alias{gs.lambda.rmla.02}
\title{Estimates or sets genomic-selection lambda values using the rmla.02 method}

\description{
Estimates or sets genomic-selection lambda values using the rmla.02 method.
}

\usage{
gs.lambda.rmla.02(alpha = 1, maxiter = 1000, precision = 0.001, hsq = 0.9,
                  out.filename = "lambda", auxfiles = TRUE,
                  data.set = "default")
}

\arguments{
   \item{alpha}{Method-specific tuning or significance parameter.}
   \item{maxiter}{Maximum number of iterations.}
   \item{precision}{Convergence tolerance.}
   \item{hsq}{Heritability parameter.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the lambda results from the RMLA.02 method when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
