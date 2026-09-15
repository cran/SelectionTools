\name{gs.mme.invcoeff}
\alias{gs.mme.invcoeff}
\title{Gs mme invcoeff}

\description{
Build or return the inverse coefficient matrix for genomic-selection mixed-model equations.
}

\usage{
gs.mme.invcoeff(out.filename = "inv.mme", auxfiles = TRUE)
}

\arguments{
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the inverse mixed-model coefficient matrix when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
