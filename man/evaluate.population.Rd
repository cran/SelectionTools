\name{evaluate.population}
\alias{evaluate.population}
\title{Evaluate genetic values for simulation populations}

\description{
Evaluates genetic values for all active individuals in one or more simulation
populations using loaded effect definitions.
}

\usage{
evaluate.population(PopNames, effName=NULL)
}

\arguments{
\item{PopNames}{Population name or population specification accepted by the
simulation routines.}
\item{effName}{Optional character string naming the effect used as the
selection index.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of calculating/storing evaluation values for the requested populations.
}

\details{
The population must have a complete genotype representation for all active
individuals. Evaluation calculates the effect-specific genetic values for all
loaded effects. Element 0 is first calculated as the weighted sum of these
values. If \code{effName} is supplied, element 0 is then replaced by the
genetic value of the named effect, preserving the existing selection
behavior.

Evaluation is population-wide. Allocation failure clears the genetic-value
category rather than leaving only some individuals evaluated.

Any stored phenotypic values for a population being re-evaluated are
discarded before the new genetic values are calculated.
}

\seealso{
\code{genotype.population}, \code{phenotype.population},
\code{get.population.gvalue}
}
