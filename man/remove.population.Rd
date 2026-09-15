\name{remove.population}
\alias{remove.population}
\title{remove.population()}

\description{
Removes one or more simulation populations and all data owned by them.
}
\usage{
remove.population(PopNames)
}

\arguments{
   \item{PopNames}{Names of the populations to be deleted}
}
\value{
Returns the same invisible implementation-level \code{list} as \code{population.remove()}. The function is a compatibility wrapper and the meaningful result is the delegated population operation.
}

\details{
Removing a population frees its chromosome representation and associated per-population storage, including cached genotype and evaluation data. Multiple population names can be supplied.
}

