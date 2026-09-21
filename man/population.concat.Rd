\name{population.concat}
\alias{population.concat}
\title{Concatenate Simulation Populations}

\description{
Moves all individuals of one population into another simulation population.
}

\usage{
population.concat(NameP1, NameP2)
}

\arguments{
\item{NameP1}{Character string naming the receiving population. If it does not
already exist, it is created.}
\item{NameP2}{Character string naming the population to concatenate and
remove.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of concatenating simulation populations.
}

\details{
Individuals already present in \code{NameP1} retain their order and remain at
the beginning of the receiving population. The active individuals from
\code{NameP2} are appended after them. Code that depends on population order
can therefore, for example, concatenate a selected fraction first and the
remaining individuals second and subsequently identify the two groups by
position.

Individuals are transferred as complete individual structures and
\code{NameP2} is removed. As in the existing implementation, cached genotype
and genetic evaluation data are discarded before concatenation. Phenotypic
values are also discarded because they depend on the genetic evaluation.
}

\seealso{
\code{population.append}, \code{population.copy}, \code{population.divide}
}
