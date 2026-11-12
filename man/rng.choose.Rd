\name{rng.choose}
\alias{rng.choose}
\title{Chooses the simulation random-number generator state}

\description{
Chooses the simulation random-number generator state.
}

\usage{
rng.choose(rngName="")
}

\arguments{
   \item{rngName}{Function argument used by this operation.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of selecting the requested random-number generator.
}
