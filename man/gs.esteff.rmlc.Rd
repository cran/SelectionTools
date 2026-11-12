\name{gs.esteff.rmlc}
\alias{gs.esteff.rmlc}
\title{gs.esteff.rmlc()}
\description{
Estimates marker effects using BLUP with constant variances for each marker.
}
\usage{
gs.esteff.rmlc(maxiter = 1000, precision = 0.0001, hsq = 0.9,
               out.filename = "eff.rmlc", auxfiles = FALSE,
               data.set ="default")
}
\arguments{
   \item{maxiter}{Maximum number of iterations used for convergence of the REML estimator}
   \item{precision}{Gives the precision as second stop value for iteration next to maxiter. If
                  the difference of the values of lambda in two subsequent iterations
                  is less than the given precision, convergence is reached and the iteration will stop}
   \item{hsq}{Function argument.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{character. Name of the data set that is used by the function}
}

\details{
The restricted maximum likelihood (REML) method is used to estimate variance components.
In a next step, lambda is calculated using the formula:

\deqn{\lambda=\frac{\sigma_e^2}{\frac{\sigma_g^2}{m}}}

where \eqn{\sigma_e^2} is the residual variance, \eqn{\sigma_g^2} is the
genetic variance and \eqn{m} is the number of markers.
}

\value{
Invisibly returns a numeric \code{matrix} containing the RMLC marker-effect estimates when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
%\note{
%
%}
%
