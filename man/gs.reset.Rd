\name{gs.reset}
\alias{gs.reset}
\title{Reset genomic-selection data and settings}

\description{
Reset genomic-selection data and settings.
}

\usage{
gs.reset()
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of resetting the genomic-selection state.
}
