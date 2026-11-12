\name{set.NoLociInit}
\alias{set.NoLociInit}
\title{Sets the initial segment capacity used by the simulation}

\description{
Sets the initial and minimum segment-storage reservation per homologue for
simulation individuals.
}

\usage{
set.NoLociInit(NoLociInit)
}

\arguments{
   \item{NoLociInit}{Positive initial segment capacity per homologue. The
   default is 10.}
}
\details{
Capacity grows geometrically when more segments are needed. Simulation
operations reserve predictable storage in advance and reuse existing capacity
on subsequent calls. This setting controls storage reservation; the stored
segments are determined by the simulation operations.

\code{population.optimize} can compact the storage of selected populations
using this reservation as the rounding unit.
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of setting the initial segment capacity.
}
