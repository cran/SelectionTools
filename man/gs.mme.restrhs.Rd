\name{gs.mme.restrhs}
\alias{gs.mme.restrhs}
\title{Gs mme restrhs}

\description{
Build or return the restricted right-hand side for genomic-selection mixed-model equations.
}

\usage{
gs.mme.restrhs(out.filename = "mme.rhs", auxfiles = FALSE)
}

\arguments{
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the restricted mixed-model right-hand side when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
