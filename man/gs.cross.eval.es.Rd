\name{gs.cross.eval.es}
\alias{gs.cross.eval.es}
\title{Evaluate crosses using expected selection response}

\description{
Evaluate crosses using expected selection response.
}

\usage{
gs.cross.eval.es(data.set = "default", alpha = 0.1, N = 0, G = 0)
}

\arguments{
   \item{data.set}{Name of the SelectionTools data set.}
   \item{alpha}{Method-specific tuning or significance parameter.}
   \item{N}{Population-size or method-specific size parameter.}
   \item{G}{Generation or method-specific generation parameter.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of evaluating crosses by expected selection response. Validation failures may return \code{NULL}.
}
