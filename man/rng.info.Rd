\name{rng.info}
\alias{rng.info}
\title{Info the simulation random-number generator state}

\description{
Info the simulation random-number generator state.
}

\usage{
rng.info()
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of requesting random-number-generator information.
}
