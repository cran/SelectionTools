\name{population.remove}
\alias{population.remove}
\title{Remove a population or population data}

\description{
Removes one or more simulation populations.
}

\usage{
population.remove(PopNames)
}

\arguments{
   \item{PopNames}{Population name or population specification.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of removing simulation populations.
}

\details{
The named populations are deleted from the simulation population registry and their owned chromosome, information, genotype, and evaluation storage is released. Other simulation-wide definitions such as the genome, linkage map, and effects are separate and are not removed by this operation.
}
