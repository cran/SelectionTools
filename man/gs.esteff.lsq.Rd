\name{gs.esteff.lsq}
\alias{gs.esteff.lsq}
\title{Estimates marker effects using the LSQ method}

\description{
Estimates marker effects using the LSQ method.
}

\usage{
gs.esteff.lsq(out.filename = "eff.lsq", auxfiles = TRUE, data.set ="default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the least-squares marker-effect estimates when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
