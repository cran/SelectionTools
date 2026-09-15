\name{population.size.get}
\alias{population.size.get}
\title{Return the size and allocated capacity of a simulation population}

\description{
Returns the active number of individuals and the number of individual slots allocated for a simulation population.
}

\usage{
population.size.get(PopName)
}

\arguments{
   \item{PopName}{Name of the population.}
}
\value{
A one-row \code{data.frame} with integer columns \code{NoInds} (active number of individuals) and \code{NoIndsAlloc} (number of individuals for which memory is allocated). Returns \code{NULL} invisibly if the population is not available.
}
