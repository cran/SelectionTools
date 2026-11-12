\name{load.internal.effmap}
\alias{load.internal.effmap}
\title{Loads an internal effect map}

\description{
Loads an internal effect map.
}

\usage{
load.internal.effmap(name, spec, weight)
}

\arguments{
   \item{name}{Function argument used by this operation.}
   \item{spec}{Function argument used by this operation.}
   \item{weight}{Weight value.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of loading the requested internal effect map.
}
