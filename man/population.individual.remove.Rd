\name{population.individual.remove}
\alias{population.individual.remove}
\title{Removes one individual from a population}

\description{
Removes one individual from a population.
}

\usage{
population.individual.remove(PopName, IndPos)
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{IndPos}{One whole-number individual position. An out-of-range position leaves the population unchanged.}
}
\value{
\code{NULL}, invisibly. A valid individual position removes that individual
while preserving the order of the remaining individuals. An out-of-range
position leaves the population unchanged.
}
