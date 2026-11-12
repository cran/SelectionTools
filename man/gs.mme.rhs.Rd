\name{gs.mme.rhs}
\alias{gs.mme.rhs}
\title{Builds or returns the right-hand side for genomic-selection mixed-model equations}

\description{
Builds or returns the right-hand side for genomic-selection mixed-model equations.
}

\usage{
gs.mme.rhs(out.filename = "mme.rhs", auxfiles = FALSE, data.set = "default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the mixed-model right-hand side when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
