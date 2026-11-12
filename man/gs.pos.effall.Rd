\name{gs.pos.effall}
\alias{gs.pos.effall}
\title{Sets or process all marker effects using the positive-effect convention}

\description{
Sets or process all marker effects using the positive-effect convention.
}

\usage{
gs.pos.effall(data.set ="default")
}

\arguments{
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of processing marker effects using the positive-effect convention. Validation failures may return \code{NULL}.
}
