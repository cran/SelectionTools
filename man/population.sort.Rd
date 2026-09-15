\name{population.sort}
\alias{population.sort}
\title{Sort a simulation population by genetic or phenotypic value}

\description{
Sorts the active individuals of a simulation population by a stored genetic
or phenotypic value.
}

\usage{
population.sort(PopName, decreasing=TRUE,
                selection.criterion="G", effect=NULL)
}

\arguments{
\item{PopName}{Character string naming the population to sort.}
\item{decreasing}{Logical. If true, larger values are placed first.}
\item{selection.criterion}{Character string. \code{"G"} selects genetic
values and is the default. \code{"P"} selects phenotypic values.}
\item{effect}{Optional character string naming an effect. If omitted, element
0 of the selected value array is used.}
}
\value{
An invisible implementation-level \code{list} returned by the compiled sorter. Its \code{retval} component is the sort status; the main effect is reordering the requested population.
}

\details{
With \code{selection.criterion="G"}, the population is sorted by
\code{GValue[0]}. With \code{selection.criterion="P"}, the corresponding
phenotypic value is used.
If \code{effect} is omitted, element 0 is the current weighted selection
index. If an effect is supplied, its effect-specific value is used.

All active individuals must have the requested value category. Sorting
reorders complete individual structures, so chromosome data, genotype data,
genetic values, phenotypic values, and information fields remain attached to
the same individual.
}

\seealso{
\code{phenotype.population}, \code{select.n.best},
\code{select.all.best}
}
