\name{set.population.gvalue}
\alias{set.population.gvalue}
\title{Sets genetic selection values for a simulation population}

\description{
Sets element 0 of the stored genetic value array for every active individual
of a simulation population.
}

\usage{
set.population.gvalue(name, gvalue)
}

\arguments{
\item{name}{Character string naming the population.}
\item{gvalue}{Finite numeric vector containing exactly one value for every
active individual.}
}
\value{
On successful assignment, an invisible implementation-level \code{list} returned by the compiled setter, including its \code{retval} status component. Invalid population, length, or non-finite input returns \code{NULL} invisibly.
}

\details{
The length of \code{gvalue} must exactly equal the active population size.
The operation is population-wide: values are set for all active individuals
or the operation is unsuccessful. Partial assignment and missing-value
padding are not used.

When effects are loaded, the R interface first genotypes and evaluates the
population so that the complete genetic value arrays exist for all effects.
The supplied values then replace element 0 only; the effect-specific genetic
values remain available. When no effects are loaded, the complete genetic
value array consists only of element 0.

A successful assignment invalidates all stored phenotypic values for the
population.
}

\seealso{
\code{get.population.gvalue}, \code{evaluate.population},
\code{phenotype.population}
}
