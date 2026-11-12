\name{gs.lambda.emstep}
\alias{gs.lambda.emstep}
\title{Estimates or sets genomic-selection lambda values using the emstep method}

\description{
Estimates or sets genomic-selection lambda values using the emstep method.
}

\usage{
gs.lambda.emstep(constvar = FALSE, out.filename = "lambda", auxfiles = TRUE,
                 data.set = "default")
}

\arguments{
   \item{constvar}{Logical flag controlling use of a constant variance.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the lambda results from the EM-step method when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
