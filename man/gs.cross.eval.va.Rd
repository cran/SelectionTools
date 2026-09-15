\name{gs.cross.eval.va}
\alias{gs.cross.eval.va}
\title{Evaluate crosses using cross variance}

\description{
Evaluate crosses using cross variance.
}

\usage{
gs.cross.eval.va(data.set = "default", pop.type = "unlinked", t = 0,
                 map.function = "Haldane")
}

\arguments{
   \item{data.set}{Name of the SelectionTools data set.}
   \item{pop.type}{Population type.}
   \item{t}{Generation or method-specific time parameter.}
   \item{map.function}{Map function to use.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of evaluating crosses by cross variance. Validation failures may return \code{NULL}.
}
