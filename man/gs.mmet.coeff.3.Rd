\name{gs.mmet.coeff.3}
\alias{gs.mmet.coeff.3}
\title{Build or solve transformed mixed-model equations (coeff.3)}

\description{
Build or solve transformed mixed-model equations (coeff.3).
}

\usage{
gs.mmet.coeff.3(out.filename = "mme.coeff", auxfiles = FALSE,
                data.set = "default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the alternative transformed mixed-model coefficient matrix when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
