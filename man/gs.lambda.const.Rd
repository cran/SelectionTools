\name{gs.lambda.const}
\alias{gs.lambda.const}
\title{Estimates or sets genomic-selection lambda values using the const method}

\description{
Estimates or sets genomic-selection lambda values using the const method.
}

\usage{
gs.lambda.const(lambda, out.filename = "lambda", auxfiles = FALSE,
                data.set = "default")
}

\arguments{
   \item{lambda}{Regularization parameter.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the lambda results from the constant-lambda calculation when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
