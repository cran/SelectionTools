\name{gs.cross.eval.gd.fct}
\alias{gs.cross.eval.gd.fct}
\title{Evaluate crosses using genetic-distance function}

\description{
Evaluate crosses using genetic-distance function.
}

\usage{
gs.cross.eval.gd.fct(dist = "rd", split = 1, data.set = "default")
}

\arguments{
   \item{dist}{Genetic-distance measure.}
   \item{split}{Split or partition control.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of evaluating factorial crosses by genetic distance. Validation failures may return \code{NULL}.
}
