\name{copy.population}
\alias{copy.population}
\title{copy.population()}

\description{
Copies all or a contiguous part of a simulation population into a destination population.
}

\usage{
copy.population(NameP1, NameP2, start=1, n=-1)
}

\arguments{
   \item{NameP1}{Name of the population to be created}
   \item{NameP2}{Name of the original population}
   \item{start}{Index of the first copied individual}
   \item{n}{Number of individuals to be copied. Default: n=-1 means all individuals}
}
\value{
Returns the same invisible implementation-level \code{list} as \code{population.copy()}. The function is a compatibility wrapper and the meaningful result is the delegated population operation.
}

\details{
Each derived cache is complete for all individuals or absent throughout a
population. A mixed source cache causes all marker-genotype, genetic-value and
phenotype caches to be cleared silently in the source and result.

Copying starts at the one-based position \code{start}. A negative \code{n} requests the remaining individuals through the end of the source; a positive value requests that many individuals subject to the source bounds. The source individuals are retained; mixed source caches are cleared as described above. Raw chromosome data are copied, and genotype or genetic-value data are copied only when the corresponding source category is complete for all active source individuals.
}


