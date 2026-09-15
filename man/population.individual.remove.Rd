\name{population.individual.remove}
\alias{population.individual.remove}
\title{Individual remove a population or population data}

\description{
Individual remove a population or population data.
}

\usage{
population.individual.remove(PopName, IndPos)
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{IndPos}{Function argument used by this operation.}
}
\value{
For a valid individual position, the invisible implementation-level \code{list} returned by the final population operation used to remove the individual. If the position is outside the population, \code{NULL}. The main effect is modification of the population.
}
