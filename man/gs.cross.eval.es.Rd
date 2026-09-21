\name{gs.cross.eval.es}
\alias{gs.cross.eval.es}
\title{Evaluate Crosses by Expected Superior Progeny Value}

\description{
Calculates the expected superior progeny value for every evaluated cross by
combining its expected mean and segregation variance.
}

\usage{
gs.cross.eval.es(data.set = "default", alpha = 0.1, N = 0, G = 0)
}

\arguments{
  \item{data.set}{Character string naming the SelectionTools data set. Cross
  means and segregation variances must already have been evaluated.}
  \item{alpha}{Selected fraction used to calculate the standardized selection
  differential when \code{N = 0} or \code{G = 0}.}
  \item{N}{Number of selected progeny for the finite-population correction.
  When both \code{N} and \code{G} are nonzero, \code{alpha} is replaced by
  \code{N/G}.}
  \item{G}{Number of evaluated progeny for the finite-population correction.
  When both \code{N} and \code{G} are nonzero, \code{alpha} is replaced by
  \code{N/G}.}
}

\details{
The expected superior progeny criterion combines the expected cross mean with
the amount of segregation expected among the progeny. For a selected fraction
\eqn{\alpha}{alpha}, the criterion described by Zhong and Jannink (2007) is

\deqn{
\mathrm{es} = \mu + i(\alpha)\,\sigma_g
}{
es = mu + i(alpha) sigma_g
}

where \eqn{\mu}{mu} is the expected cross mean and
\eqn{\sigma_g}{sigma_g} is the square root of the segregation variance. The
standardized selection differential is

\deqn{
i(\alpha) =
\frac{g_{0,1}\left[G^{-1}_{0,1}(1-\alpha)\right]}{\alpha}
}{
i(alpha) = g_0,1[G^-1_0,1(1-alpha)] / alpha
}

where \eqn{g_{0,1}}{g_0,1} is the probability density function of the standard
normal distribution and \eqn{G^{-1}_{0,1}}{G^-1_0,1} is the inverse of its
cumulative distribution function.

For a small finite progeny population, when both \code{N} and \code{G} are
nonzero, the function sets \eqn{\alpha=N/G}{alpha=N/G} and applies the
correction of Burrows (1972):

\deqn{
\mathrm{es} = \mu + i(N,G)\,\sigma_g
}{
es = mu + i(N,G) sigma_g
}

with

\deqn{
i(N,G) = i(\alpha) -
\frac{G-N}{2N(G+1)i(\alpha)}
}{
i(N,G) = i(alpha) - (G-N) / [2 N (G+1) i(alpha)]
}

The calculated value is stored as the cross criterion \code{es}. Results for
all evaluated crosses can subsequently be obtained with \code{gs.cross.info}.
}

\references{
Zhong S, Jannink JL (2007) Using quantitative trait loci results to
discriminate among crosses on the basis of their progeny mean and variance.
Genetics 177:567--576.

Burrows RM (1972) Expected selection differentials for directional selection.
Biometrics 28:1091--1100.
}

\value{
On a successful run, an invisible implementation-level \code{list} returned
by the final timing/information helper; the meaningful result is the side
effect of evaluating all crosses by expected superior progeny value.
Validation failures may return \code{NULL}.
}

\seealso{
\code{gs.cross.eval.mu}, \code{gs.cross.eval.va}, \code{gs.cross.info}
}
