\name{gs.vc.rrblup}
\alias{gs.vc.rrblup}
\title{Estimates variance components for ridge-regression BLUP}

\description{
Estimates variance components for ridge-regression BLUP.
}

\usage{
gs.vc.rrblup(maxiter = 1000, precision = 1e-4, hsq = 0.8,
             out.filename = "lambda", auxfiles = TRUE, data.set = "default")
}

\arguments{
   \item{maxiter}{Main optimizer iteration budget. The shared engine may
   perform additional convergence polishing beyond this budget.}
   \item{precision}{Projected-gradient convergence tolerance in the shared
   mixed-model engine's scaled covariance parameters.}
   \item{hsq}{Heritability parameter.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the variance-component results from RR-BLUP when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
