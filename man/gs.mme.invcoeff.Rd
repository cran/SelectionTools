\name{gs.mme.invcoeff}
\alias{gs.mme.invcoeff}
\title{Gs mme invcoeff}

\description{
Builds or returns the inverse coefficient matrix for genomic-selection mixed-model equations.
}

\usage{
gs.mme.invcoeff(out.filename = "inv.mme", auxfiles = TRUE)
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the inverse mixed-model coefficient matrix when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
