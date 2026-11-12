\name{population.name.swap}
\alias{population.name.swap}
\title{Swap the names of two simulation populations}

\description{
Swaps the supported names of two simulation populations.
}

\usage{
population.name.swap(NameP1, NameP2)
}

\arguments{
   \item{NameP1}{Name of the first parental population.}
   \item{NameP2}{Name of the second parental population.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of swapping simulation-population names.
}
