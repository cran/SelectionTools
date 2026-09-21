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
   \item{PopName}{Character string containing one population name or several
population names separated by blanks.}
}

\details{
The active size is the number of individuals currently belonging to a
population. The allocated capacity is internal storage available to the
population and can differ from the active size after resizing or other
population operations.

Several population names can be supplied in one character string, separated
by blanks. In that case, \code{NoInds} and \code{NoIndsAlloc} are summed
over all specified populations.
}
\value{
Returns the same one-row \code{data.frame} as \code{population.size.get()},
with columns \code{NoInds} and \code{NoIndsAlloc}. For several population
names these columns contain the corresponding sums. Returns \code{NULL} if
a requested population is unavailable.
}

