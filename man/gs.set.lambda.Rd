\name{gs.set.lambda}
\alias{gs.set.lambda}
\title{Set lambda used by the genomic-selection routines}

\description{
Set lambda used by the genomic-selection routines.
}

\usage{
gs.set.lambda(lambda, out.filename = "lambda", auxfiles = FALSE,
              data.set = "default")
}

\arguments{
   \item{lambda}{Regularization parameter.}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the lambda values produced while setting lambda when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
