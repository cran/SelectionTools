\name{set.NoLociInit}
\alias{set.NoLociInit}
\title{Set the initial number of loci used by the simulation}

\description{
Set the initial number of loci used by the simulation.
}

\usage{
set.NoLociInit(NoLociInit)
}

\arguments{
   \item{NoLociInit}{Initial number of loci.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of setting the initial number of loci.
}
