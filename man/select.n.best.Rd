\name{select.n.best}
\alias{select.n.best}
\title{Selects a fixed number of best individuals}

\description{
Selects a fixed number of individuals from a simulation population according
to genetic or phenotypic values.
}

\usage{
select.n.best(newPop, oldPop, effect=NULL, n=1,
              decreasing=TRUE, selection.criterion="G")
}

\arguments{
\item{newPop}{Character string naming the selected population.}
\item{oldPop}{Character string naming the source population.}
\item{effect}{Optional character string naming the effect used for selection.}
\item{n}{Number of individuals requested.}
\item{decreasing}{Logical. If true, larger values are preferred.}
\item{selection.criterion}{Character string. \code{"G"}, the default,
selects on genetic values; \code{"P"} selects on phenotypic values.}
}
\value{
On successful selection, returns the invisible implementation-level \code{list} produced by \code{population.divide()}; if sorting fails the sort-status list is returned invisibly, and an empty population returns \code{NULL} invisibly. The main effect is creation of \code{newPop} containing the selected individuals.
}

\details{
Each derived cache is complete for all individuals or absent throughout a
population. A mixed source cache causes all marker-genotype, genetic-value and
phenotype caches to be cleared silently in the source and result.

With \code{selection.criterion="G"}, the source population is genotyped and
evaluated before selection. With \code{effect=NULL}, selection is based on
the weighted genetic index.
With a named \code{effect}, selection is based on that effect.

With \code{selection.criterion="P"}, stored phenotypic values are used and
the population is not re-evaluated, because genetic re-evaluation invalidates
phenotypic values. If \code{effect=NULL}, the weighted phenotypic index is
used; with a named effect, its effect-specific phenotype is used.
}

\seealso{
\code{phenotype.population}, \code{population.sort},
\code{select.all.best}
}
