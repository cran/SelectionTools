\name{gs.mme.restrhs}
\alias{gs.mme.restrhs}
\title{Gs mme restrhs}

\description{
Builds or returns the restricted right-hand side for genomic-selection mixed-model equations.
}

\usage{
gs.mme.restrhs(out.filename = "mme.rhs", auxfiles = FALSE)
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the restricted mixed-model right-hand side when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
