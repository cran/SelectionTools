\name{gs.mme.restcoeff}
\alias{gs.mme.restcoeff}
\title{Gs mme restcoeff}

\description{
Builds or returns the penalized coefficient matrix for genomic-selection mixed-model equations.
}

\usage{
gs.mme.restcoeff(out.filename = "mme.coeff", auxfiles = FALSE)
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the penalized mixed-model coefficient matrix when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
