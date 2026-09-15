\name{population.optimize}
\alias{population.optimize}
\title{Optimize storage for simulation populations}

\description{
Calls the compiled population-optimization routine for the requested simulation populations.
}

\usage{
population.optimize(PopNames)
}

\arguments{
   \item{PopNames}{Population name or population specification.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of optimizing the requested simulation populations.
}
