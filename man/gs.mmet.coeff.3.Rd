\name{gs.mmet.coeff.3}
\alias{gs.mmet.coeff.3}
\title{Builds or solves transformed mixed-model equations (coeff.3)}

\description{
Builds or solves transformed mixed-model equations (coeff.3).
}

\usage{
gs.mmet.coeff.3(out.filename = "mme.coeff", auxfiles = FALSE,
                data.set = "default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the alternative transformed mixed-model coefficient matrix when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
