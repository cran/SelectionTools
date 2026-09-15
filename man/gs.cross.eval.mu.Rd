\name{gs.cross.eval.mu}
\alias{gs.cross.eval.mu}
\title{Evaluate crosses using cross mean}

\description{
Evaluate crosses using cross mean.
}

\usage{
gs.cross.eval.mu(data.set = "default")
}

\arguments{
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of evaluating crosses by expected cross mean. Validation failures may return \code{NULL}.
}
