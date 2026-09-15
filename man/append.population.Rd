\name{append.population}
\alias{append.population}
\title{append.population()}

\description{
Appends a copy of one population to another.
}

\usage{
append.population(NameP1, NameP2)
}

\arguments{
   \item{NameP1}{Population which takes up the copied individuals}
   \item{NameP2}{The population to be appended}
}
\value{
Returns the same invisible implementation-level \code{list} as \code{population.append()}. The function is a compatibility wrapper and the meaningful result is the delegated population operation.
}

\details{
The individuals of \code{NameP2} are copied to the end of \code{NameP1}; \code{NameP2} is not modified. Raw chromosome data are copied for every appended individual. Cached genotype and genetic-value categories are retained only when they can remain complete and internally consistent for all active individuals of the destination; otherwise the corresponding derived category is cleared.
}

