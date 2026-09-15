\name{population.divide}
\alias{population.divide}
\title{Divide a simulation population}

\description{
Moves the first requested individuals of a source population into a new
population.
}

\usage{
population.divide(NameP1, NameP2, NoI=-1)
}

\arguments{
\item{NameP1}{Character string naming the destination population.}
\item{NameP2}{Character string naming the source population.}
\item{NoI}{Number of individuals to move.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of dividing a simulation population.
}

\details{
The operation moves complete individual structures. Therefore chromosome
data, information fields, genotype data, genetic values, and phenotypic
values remain attached to the same individuals. No value category is
recalculated merely because the population is divided.
}
