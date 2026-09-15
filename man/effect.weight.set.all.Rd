\name{effect.weight.set.all}
\alias{effect.weight.set.all}
\title{Set effect weights for all effects}

\description{
Set effect weights for all effects.
}

\usage{
effect.weight.set.all(weight)
}

\arguments{
   \item{weight}{Weight value.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of setting all effect weights.
}
