\name{gs.set.lambda}
\alias{gs.set.lambda}
\title{Sets lambda used by the genomic-selection routines}

\description{
Sets lambda used by the genomic-selection routines.
}

\usage{
gs.set.lambda(lambda, out.filename = "lambda", auxfiles = FALSE,
              data.set = "default")
}

\arguments{
   \item{lambda}{Regularization parameter.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the lambda values produced while setting lambda when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
