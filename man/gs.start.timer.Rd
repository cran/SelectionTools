\name{gs.start.timer}
\alias{gs.start.timer}
\title{Starts the genomic-selection timing helper}

\description{
Starts the genomic-selection timing helper.
}

\usage{
gs.start.timer(depth=1)
}

\arguments{
   \item{depth}{Timer depth, a finite whole number from 1 through 5.}
}
\value{
An invisible implementation-level \code{list} returned by the SelectionTools timing helper. The function is called for the side effect of starting timing.
}
