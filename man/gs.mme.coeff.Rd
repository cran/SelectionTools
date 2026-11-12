\name{gs.mme.coeff}
\alias{gs.mme.coeff}
\title{Gs mme coeff}

\description{
Builds or returns the coefficient matrix for genomic-selection mixed-model equations.
}

\usage{
gs.mme.coeff(out.filename = "mme.coeff", auxfiles = FALSE,
             data.set = "default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the mixed-model coefficient matrix when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
