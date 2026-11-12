\name{gs.cross.eval.mu}
\alias{gs.cross.eval.mu}
\title{Evaluates Crosses by Expected Progeny Mean}

\description{
Calculates the expected mean performance of the progeny for every possible
pair of parents in a SelectionTools data set from previously estimated
marker effects.
}

\usage{
gs.cross.eval.mu(data.set = "default")
}

\arguments{
  \item{data.set}{Character string naming the SelectionTools data set. The
  data set must contain estimated marker effects.}
}

\details{
For a cross between parents P1 and P2, the expected progeny mean is
calculated as

\deqn{
\mu = \beta_0 + \sum_m (\beta_{mi} + \beta_{mj})/2
}{
mu = beta0 + sum_m (beta_mi + beta_mj)/2
}

where \eqn{\beta_0}{beta0} is the intercept,
\eqn{\beta_{mi}}{beta_mi} is the effect of allele \eqn{i} carried by P1 at
marker \eqn{m}, and \eqn{\beta_{mj}}{beta_mj} is the effect of allele
\eqn{j} carried by P2 at marker \eqn{m}.

The calculated value is stored as the cross criterion \code{mu}. Results for
all evaluated crosses can subsequently be obtained with \code{gs.cross.info}.
}

\value{
On a successful run, an invisible implementation-level \code{list} returned
by the final timing/information helper; the meaningful result is the side
effect of evaluating all crosses by expected progeny mean. Validation
failures may return \code{NULL}.
}

\seealso{
\code{gs.cross.eval.mi}, \code{gs.cross.eval.ma},
\code{gs.cross.eval.va}, \code{gs.cross.eval.es}, \code{gs.cross.info}
}
