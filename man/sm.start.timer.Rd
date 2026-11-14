\name{sm.start.timer}
\alias{sm.start.timer}
\title{Starts the simulation timing helper}

\description{
Starts the simulation timing helper.
}

\usage{
sm.start.timer(depth=1)
}

\arguments{
   \item{depth}{Timer depth, a finite whole number from 1 through 5.}
}
\value{
An invisible implementation-level \code{list} returned by the SelectionTools timing helper. The function is called for the side effect of starting timing.
}
