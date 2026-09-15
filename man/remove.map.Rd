\name{remove.map}
\alias{remove.map}
\title{Remove the linkage map from the simulation}

\description{
Removes the current linkage map and map-dependent cached simulation data.
}

\usage{
remove.map()
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of removing the current linkage map.
}

\details{
Removing the linkage map clears cached marker genotypes, genetic values,
phenotypic values, and loaded effect definitions because these data depend on
the current map/effect configuration. Population chromosome-segment data are
retained according to the simulation backend.
}
