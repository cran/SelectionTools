\name{get.population.size}
\alias{get.population.size}
\title{get.population.size()}

\description{
Returns the active size and allocated capacity of a simulation population.
}

\usage{
get.population.size(PopName)
}

\arguments{
   \item{PopName}{Name of the population}
}

\details{
The active size is the number of individuals currently belonging to the population. The allocated capacity is internal storage available to the population and can differ from the active size after resizing or other population operations.
}
\value{
Returns the same one-row \code{data.frame} as \code{population.size.get()}, with columns \code{NoInds} and \code{NoIndsAlloc}, or \code{NULL} if the population is unavailable.
}

