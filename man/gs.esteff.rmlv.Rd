% copyright (C) 2001 Frisch/Maurer
%
\name{gs.esteff.rmlv}
\alias{gs.esteff.rmlv}
\title{gs.esteff.rmlv()}
\description{
Estimates marker effects using BLUP with variable variances for the markers.
}
\usage{
gs.esteff.rmlv(maxiter = 1000, precision = 0.001, hsq = 0.9,
               out.filename = "eff.rmlv", auxfiles = FALSE,
               data.set ="default")
}
\arguments{
   \item{maxiter}{Maximum number of iterations used for convergence of the REML estimator.}
   \item{precision}{Gives the precision as second stop value for iteration next to maxiter. If
                  the difference of the values of lambda in two subsequent iterations
                  is less than the given precision, convergence is reached and the iteration will stop.}
   \item{hsq}{Function argument.}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{character. Name of the data set that is used by the function.}
}

\details{
The restricted maximum likelihood (REML) method is used to estimate variance components.

}

\value{
Invisibly returns a numeric \code{matrix} containing the RMLV marker-effect estimates when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}

%\note{
%
%}
%
