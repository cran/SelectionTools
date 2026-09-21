\name{gs.esteff.rr}
\alias{gs.esteff.rr}
\title{Estimate Marker Effects by Ridge Regression}

\description{
Estimates genome-wide marker effects by ridge-regression methods, including
RR-BLUP.
}

\usage{
gs.esteff.rr(method = "BLUP", maxiter = 1000, precision = 0.0001, hsq = 0.80,
             alpha = 1, out.filename = "eff.rr", auxfiles = FALSE,
             data.set = "default")
}

\arguments{
   \item{method}{Character string selecting the ridge-regression method. The
   default \code{"BLUP"} performs RR-BLUP.}
   \item{maxiter}{Maximum number of iterations used for variance-component
   estimation where applicable. For \code{method="BLUP"}, setting
   \code{maxiter=0} skips variance-component iteration and uses the variance
   components implied by \code{hsq} directly.}
   \item{precision}{Convergence tolerance for iterative variance-component
   estimation.}
   \item{hsq}{Heritability used to construct starting variance components. For
   \code{method="BLUP"} with \code{maxiter>0}, these values initialize the
   variance-component iteration. With \code{maxiter=0}, they are used without
   further variance-component estimation.}
   \item{alpha}{Method-specific tuning or significance parameter.}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set containing marker and
   performance data.}
}

\details{
With \code{method="BLUP"}, SelectionTools fits the standard ridge-regression
BLUP model. The implementation derives starting genetic and residual variance
components from the observed performance variance and \code{hsq}. If
\code{maxiter>0}, these starting values are refined by the mixed-model engine
until the requested precision is reached or \code{maxiter} iterations have
been performed.

Setting \code{maxiter=0} omits this variance-component iteration. In that
case the \code{hsq}-derived variance components are retained directly. This
provides the computationally faster fixed-heritability form used in the
SelectionTools simulation vignette.

After effect estimation, the marker effects are stored in \code{data.set} and
can be used for genomic prediction, cross evaluation, or transferred to the
simulation backend. The model fit for the training population can be inspected
with \code{gs.plot.model.fit}.
}

\value{
Invisibly returns a numeric \code{matrix} containing the estimated marker effects from ridge regression when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}. The estimated effects are also stored in \code{data.set}.
}

\references{
Whittaker JC, Thompson R, Denham MC (2000) Marker-assisted selection using
ridge regression. Genetics Research 75:249--252.

Meuwissen THE, Hayes BJ, Goddard ME (2001) Prediction of total genetic value
using genome-wide dense marker maps. Genetics 157:1819--1829.

Habier D, Fernando RL, Dekkers JCM (2007) The impact of genetic relationship
information on genome-assisted breeding values. Genetics 177:2389--2397.

Hofheinz N, Frisch M (2014) Heteroscedastic ridge regression approaches for
genome-wide prediction with a focus on computational efficiency and accurate
effect estimation. G3 4:539--546.
}

\seealso{
\code{gs.plot.model.fit}, \code{gs.predict.genotypes}, \code{gs.return.effects}
}
