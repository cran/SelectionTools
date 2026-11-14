\name{optimize.population}
\alias{optimize.population}
\title{Optimizes storage for simulation populations}

\description{
Compacts segment storage for the requested simulation populations while
preserving all stored segments.
}

\usage{
optimize.population(PopNames)
}

\arguments{
   \item{PopNames}{Population name or population specification.}
}
\value{
Returns the same invisible implementation-level \code{list} as \code{population.optimize()}. The function is a compatibility wrapper and the meaningful result is the delegated population operation.
}

\details{
Segment storage is compacted for the actual individuals in each requested
population. Unused individual slots are released and the individual array is
trimmed to the population size. Chromosomes, individual identities and complete
derived caches are preserved. Segment capacity follows the reservation set by
\code{set.NoLociInit}.

Completely free backing blocks are released from the shared simulation memory
pool after successful optimization. Blocks containing any live allocation are
retained. Other populations' live allocations are not moved or released.
Released memory need not immediately reduce the process's resident memory on
every system. Later simulation growth can allocate storage again.
}
