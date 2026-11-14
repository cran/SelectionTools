\name{divide.population}
\alias{divide.population}
\title{divide.population()}

\description{
Splits individuals from one simulation population into a second population.
}

\usage{
divide.population(NameP1, NameP2, NoI=-1)
}

\arguments{
   \item{NameP1}{Name of the population to be created}
   \item{NameP2}{Name of the population to be divided}
   \item{NoI}{Number of individuals to be separated from population \emph{NameP2}}
}
\value{
Returns the same invisible implementation-level \code{list} as \code{population.divide()}. The function is a compatibility wrapper and the meaningful result is the delegated population operation.
}

\details{
Each derived cache is complete for all individuals or absent throughout a
population. A mixed source cache causes all marker-genotype, genetic-value and
phenotype caches to be cleared silently in the source and result.

  The first requested individuals are
  separated from \code{NameP2} and placed in \code{NameP1}; the
  remaining individuals stay in \code{NameP2}. Complete individual
  structures are moved, so chromosome segments, information fields,
  genotype data, and evaluations that belong to an individual move
  together.
}

