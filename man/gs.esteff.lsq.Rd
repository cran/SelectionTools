\name{gs.esteff.lsq}
\alias{gs.esteff.lsq}
\title{Estimate marker effects using the LSQ method}

\description{
Estimate marker effects using the LSQ method.
}

\usage{
gs.esteff.lsq(out.filename = "eff.lsq", auxfiles = TRUE, data.set ="default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the least-squares marker-effect estimates when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
