\name{resize.population}
\alias{resize.population}
\title{Resize a population}

\description{
Changes the allocated size of a simulation population.
}

\usage{
resize.population(PopName, newSize)
}

\arguments{
\item{PopName}{Character string naming the population to resize.}
   \item{newSize}{Requested new population allocation.}
}
\value{
Returns the same invisible implementation-level \code{list} as \code{population.resize()}. The function is a compatibility wrapper and the meaningful result is the delegated population operation.
}

\details{
A successful change of allocation invalidates cached marker genotypes and evaluated genetic values because these derived categories must remain complete and consistent for the active population. If no actual resize occurs, the cached derived data are left unchanged. If an attempted growth fails and the original allocation remains unchanged, the original population state is preserved.
}

