\name{gs.lambda.const}
\alias{gs.lambda.const}
\title{Estimate or set genomic-selection lambda values using the const method}

\description{
Estimate or set genomic-selection lambda values using the const method.
}

\usage{
gs.lambda.const(lambda, out.filename = "lambda", auxfiles = FALSE,
                data.set = "default")
}

\arguments{
   \item{lambda}{Regularization parameter.}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the lambda results from the constant-lambda calculation when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
