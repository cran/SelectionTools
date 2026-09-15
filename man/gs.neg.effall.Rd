\name{gs.neg.effall}
\alias{gs.neg.effall}
\title{Set or process all marker effects using the negative-effect convention}

\description{
Set or process all marker effects using the negative-effect convention.
}

\usage{
gs.neg.effall(data.set ="default")
}

\arguments{
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of processing marker effects using the negative-effect convention. Validation failures may return \code{NULL}.
}
