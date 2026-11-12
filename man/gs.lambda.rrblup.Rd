\name{gs.lambda.rrblup}
\alias{gs.lambda.rrblup}
\title{Estimates or sets genomic-selection lambda values using the rrblup method}

\description{
Estimates or sets genomic-selection lambda values using the rrblup method.
}

\usage{
gs.lambda.rrblup(maxiter = 1000, precision = 1e-4, hsq = 0.8,
                 out.filename = "lambda", auxfiles = TRUE,
                 data.set = "default")
}

\arguments{
   \item{maxiter}{Maximum number of iterations.}
   \item{precision}{Convergence tolerance.}
   \item{hsq}{Heritability parameter.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the lambda results from RR-BLUP when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
