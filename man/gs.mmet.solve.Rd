\name{gs.mmet.solve}
\alias{gs.mmet.solve}
\title{Builds or solves transformed mixed-model equations (solve)}

\description{
Builds or solves transformed mixed-model equations (solve).
}

\usage{
gs.mmet.solve(out.filename = "mme.solution", auxfiles = FALSE,
              data.set = "default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the transformed mixed-model solution when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
