\name{gs.cross.eval.mi}
\alias{gs.cross.eval.mi}
\title{Evaluate crosses using marker information}

\description{
Evaluate crosses using marker information.
}

\usage{
gs.cross.eval.mi(data.set = "default")
}

\arguments{
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of evaluating crosses by marker information. Validation failures may return \code{NULL}.
}
