\name{gs.cross.eval.ma}
\alias{gs.cross.eval.ma}
\title{Evaluate Crosses by Maximum Marker-Effect Value}

\description{
Calculates the maximum marker-effect based performance that can be assembled
from the parental alleles for every possible pair of parents in a
SelectionTools data set.
}

\usage{
gs.cross.eval.ma(data.set = "default")
}

\arguments{
  \item{data.set}{Character string naming the SelectionTools data set. The
  data set must contain estimated marker effects.}
}

\details{
Using the notation of \code{gs.cross.eval.mu}, the maximum value for a line
derived from a cross is

\deqn{
\mathrm{ma} = \beta_0 + \sum_m \max(\beta_{mi},\beta_{mj})
}{
ma = beta0 + sum_m max(beta_mi,beta_mj)
}

where \eqn{\beta_0}{beta0} is the intercept and
\eqn{\beta_{mi}}{beta_mi} and \eqn{\beta_{mj}}{beta_mj} are the marker-effect
contributions of the alleles carried by P1 and P2 at marker \eqn{m}.

The calculated value is stored as the cross criterion \code{ma}. It is a
marker-effect based upper limit obtained by combining, marker by marker, the
larger of the two parental contributions. Results for all evaluated crosses
can subsequently be obtained with \code{gs.cross.info}.
}

\value{
On a successful run, an invisible implementation-level \code{list} returned
by the final timing/information helper; the meaningful result is the side
effect of evaluating all crosses by the maximum marker-effect value.
Validation failures may return \code{NULL}.
}

\seealso{
\code{gs.cross.eval.mu}, \code{gs.cross.eval.mi}, \code{gs.cross.info}
}
