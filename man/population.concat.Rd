\name{population.concat}
\alias{population.concat}
\title{Concatenate simulation populations}

\description{
Moves all individuals of one population into another simulation population.
}

\usage{
population.concat(NameP1, NameP2)
}

\arguments{
\item{NameP1}{Character string naming the receiving population.}
\item{NameP2}{Character string naming the population to concatenate and
remove.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of concatenating simulation populations.
}

\details{
Individuals are transferred as complete individual structures and
\code{NameP2} is removed. As in the existing implementation, cached genotype
and genetic evaluation data are discarded before concatenation. Phenotypic
values are also discarded because they depend on the genetic evaluation.
}
