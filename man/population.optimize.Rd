\name{population.optimize}
\alias{population.optimize}
\title{Optimize storage for simulation populations}

\description{
Compacts segment storage for the requested simulation populations while
preserving all stored segments.
}

\usage{
population.optimize(PopNames)
}

\arguments{
   \item{PopNames}{Population name or population specification.}
}
\details{
Ordinary simulation operations retain allocated capacity for reuse. This
function explicitly compacts each individual's segment buffer, keeping the
segments of every allocated individual slot, including temporarily inactive
slots. Segment capacity is rounded using the reservation set by
\code{set.NoLociInit}.

Released storage becomes available to other simulation populations in the
package's compiled memory pool. The pool retains its backing blocks for later
operations, so compaction need not reduce the process's resident memory.
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of optimizing the requested simulation populations.
}
