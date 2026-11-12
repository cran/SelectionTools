\name{remove.evaluate.population}
\alias{remove.evaluate.population}
\title{Removes genetic and phenotypic evaluation values}

\description{
Clears stored genetic and phenotypic values for one or more simulation
populations.
}

\usage{
remove.evaluate.population(PopNames)
}

\arguments{
\item{PopNames}{Population name or population specification.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of removing stored evaluation values.
}

\details{
Population membership, chromosome segments, and cached marker genotypes are
retained. Both genetic values and phenotypic values are discarded. Genetic
values must be evaluated or set again before genetic selection, and
phenotypes must be generated again before phenotypic selection.
}
