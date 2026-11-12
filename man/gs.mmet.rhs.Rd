\name{gs.mmet.rhs}
\alias{gs.mmet.rhs}
\title{Builds or solves transformed mixed-model equations (rhs)}

\description{
Builds or solves transformed mixed-model equations (rhs).
}

\usage{
gs.mmet.rhs(out.filename = "mme.rhs", auxfiles = FALSE, data.set = "default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the transformed mixed-model right-hand side when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
