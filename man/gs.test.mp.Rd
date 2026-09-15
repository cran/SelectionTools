\name{gs.test.mp}
\alias{gs.test.mp}
\title{Test marker or model parameters used by genomic-selection routines}

\description{
Test marker or model parameters used by genomic-selection routines.
}

\usage{
gs.test.mp(n = 0)
}

\arguments{
   \item{n}{Number of items or repetitions.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of running the compiled marker/model-parameter test.
}
