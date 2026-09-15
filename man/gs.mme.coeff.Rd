\name{gs.mme.coeff}
\alias{gs.mme.coeff}
\title{Gs mme coeff}

\description{
Build or return the coefficient matrix for genomic-selection mixed-model equations.
}

\usage{
gs.mme.coeff(out.filename = "mme.coeff", auxfiles = FALSE,
             data.set = "default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the mixed-model coefficient matrix when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
