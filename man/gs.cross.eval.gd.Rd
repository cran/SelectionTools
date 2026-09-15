\name{gs.cross.eval.gd}
\alias{gs.cross.eval.gd}
\title{Evaluate crosses using genetic distance}

\description{
Evaluate crosses using genetic distance.
}

\usage{
gs.cross.eval.gd(dist = "rd", data.set = "default")
}

\arguments{
   \item{dist}{Genetic-distance measure.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of evaluating crosses by genetic distance. Validation failures may return \code{NULL}.
}
