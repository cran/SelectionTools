\name{gs.stop.timer}
\alias{gs.stop.timer}
\title{Stops the genomic-selection timing helper}

\description{
Stops the genomic-selection timing helper.
}

\usage{
gs.stop.timer(depth=1, info.level=0)
}

\arguments{
   \item{depth}{Timer depth, a finite whole number from 1 through 5.}
   \item{info.level}{Function argument used by this operation.}
}
\value{
An invisible implementation-level \code{list} returned by the SelectionTools information helper after reporting elapsed time. The function is called for the timing/reporting side effect.
}
