\name{init.population}
\alias{init.population}
\title{init.population()}

\description{
Creates a simulation population from the internal chromosome-segment data-frame representation.
}

\usage{
init.population(name, data, delete=TRUE)
}

\arguments{
\item{name}{Character string naming the simulation population to initialize.}
   \item{data}{Data frame with five columns representing individual, chromosome, homologue, segment position, and allele.}
   \item{delete}{Logical. If true, remove an existing population with the same name before initialization.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of initializing the simulation population from the supplied data.
}

\details{
The input data frame must have exactly five columns named or interpretable as individual number, chromosome number, homologue number, position, and allele. Each row describes the start of a chromosome segment carrying the specified allele. The data must be ordered consistently with the simulation representation.

The number of input rows must be positive. Chromosome and homologue identifiers must be compatible with the currently defined simulation genome. The routine constructs the raw chromosome representation; marker genotypes and evaluated genetic values are derived data and are not implied merely by initialization.

When \code{delete=TRUE}, an existing population of the same name is removed before initialization. If creation fails while chromosome segments are being inserted, the incomplete new population is removed rather than retained.
}
