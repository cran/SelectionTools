\name{gs.mme.solve}
\alias{gs.mme.solve}
\title{Gs mme solve}

\description{
Builds or returns the mixed-model equation solution for genomic-selection mixed-model equations.
}

\usage{
gs.mme.solve(out.filename = "mme.solution", auxfiles = FALSE,
             data.set = "default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the mixed-model equation solution when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
