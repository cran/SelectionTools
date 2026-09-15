\name{sample.population}
\alias{sample.population}
\title{sample.population()}

\description{
Creates a new population by sampling individuals from an existing simulation population.
}

\usage{
sample.population(NameP1, NameP2, NoI=-1, rep=1)
}

\arguments{
   \item{NameP1}{Name of the sample population}
   \item{NameP2}{Name of the initial population}
   \item{NoI}{non-negative integer giving the number of individuals to choose}
   \item{rep}{Should sampling be with replacement?}
}
\value{
Returns the same invisible implementation-level \code{list} as \code{population.sample()}. The function is a compatibility wrapper and the meaningful result is the delegated population operation.
}

\details{
Sampling does not modify the source population. The requested sample is copied into the destination. Sampling can be with or without replacement according to \code{rep}. Raw chromosome data are always the authoritative copied representation; cached genotype and genetic-value categories are transferred only when complete in the source and successfully constructed for the entire destination.
}

