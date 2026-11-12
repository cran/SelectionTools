\name{population.remove.all}
\alias{population.remove.all}
\title{Removes all simulation populations}

\description{
Removes all currently supported simulation populations.
}

\usage{
population.remove.all()
}
\value{
The function returns \code{NULL} invisibly after the population-removal request.
}

\details{
The function obtains the current population list and removes every listed population. This clears population-owned data but does not by itself redefine genome parameters, linkage maps, or effect maps.

Population storage is returned to the package's compiled memory pool. Backing
blocks remain available for subsequent simulations, as described for
\code{population.remove}.
}
