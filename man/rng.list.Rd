\name{rng.list}
\alias{rng.list}
\title{List the simulation random-number generator state}

\description{
Lists the available simulation random-number generators and indicate the current
choice. Output follows the package information-level setting and can therefore
be suppressed through that setting.
}

\usage{
rng.list()
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of requesting the available random-number generators.
}
