\name{gs.mmet.solve}
\alias{gs.mmet.solve}
\title{Build or solve transformed mixed-model equations (solve)}

\description{
Build or solve transformed mixed-model equations (solve).
}

\usage{
gs.mmet.solve(out.filename = "mme.solution", auxfiles = FALSE,
              data.set = "default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the transformed mixed-model solution when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
