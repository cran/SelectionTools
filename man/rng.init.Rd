\name{rng.init}
\alias{rng.init}
\title{Init the simulation random-number generator state}

\description{
Init the simulation random-number generator state.
}

\usage{
rng.init()
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of initializing the random-number generator.
}
