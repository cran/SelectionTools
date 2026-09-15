\name{population.swap.name}
\alias{population.swap.name}
\title{Swap simulation population names}

\description{
Compatibility wrapper for swapping the names of two simulation populations.
}

\usage{
population.swap.name(NameP1, NameP2)
}

\arguments{
   \item{NameP1}{Name of the first parental population.}
   \item{NameP2}{Name of the second parental population.}
}
\value{
Returns the same invisible implementation-level \code{list} as \code{population.name.swap()}. The function is a compatibility wrapper and the meaningful result is the delegated population operation.
}
