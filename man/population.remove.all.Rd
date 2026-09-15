\name{population.remove.all}
\alias{population.remove.all}
\title{Remove all a population or population data}

\description{
Removes all currently registered simulation populations.
}

\usage{
population.remove.all()
}
\value{
An invisible implementation-level \code{list} returned by the final information-message call after all populations have been removed. The meaningful result is the side effect of removing all populations.
}

\details{
The function obtains the current population list and removes every listed population. This clears population-owned data but does not by itself redefine genome parameters, linkage maps, or effect maps.
}
