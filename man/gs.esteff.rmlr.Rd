\name{gs.esteff.rmlr}
\alias{gs.esteff.rmlr}
\title{Estimates marker effects using the RMLR method}

\description{
Estimates marker effects using the RMLR method.
}

\usage{
gs.esteff.rmlr(maxiter = 1000, precision = 0.001, hsq = 0.9,
               out.filename = "eff.rmlr", auxfiles = TRUE,
               data.set ="default")
}

\arguments{
   \item{maxiter}{Maximum number of iterations.}
   \item{precision}{Convergence tolerance.}
   \item{hsq}{Heritability parameter.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the RMLR marker-effect estimates when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
