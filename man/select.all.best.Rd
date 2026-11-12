\name{select.all.best}
\alias{select.all.best}
\title{Selects all individuals in the best value classes}

\description{
Selects individuals belonging to the best distinct genetic or phenotypic
value classes of a simulation population.
}

\usage{
select.all.best(newPop, oldPop, effect=NULL, x=1,
                decreasing=TRUE, selection.criterion="G")
}

\arguments{
\item{newPop}{Character string naming the selected population.}
\item{oldPop}{Character string naming the source population.}
\item{effect}{Optional character string naming the effect used for selection.}
\item{x}{Number of distinct value classes to retain.}
\item{decreasing}{Logical. If true, larger values are preferred.}
\item{selection.criterion}{Character string. \code{"G"}, the default,
selects on genetic values; \code{"P"} selects on phenotypic values.}
}
\value{
On successful selection, invisibly returns a one-row \code{data.frame} with columns \code{n}, \code{minscore}, and \code{maxscore}, giving the number selected and the score limits of the selected class(es). Returns \code{NULL} for an empty/unusable selection, or the invisible status \code{list} from \code{population.sort()} if sorting fails.
}

\details{
With \code{selection.criterion="G"}, the source population is genotyped and
evaluated before sorting. Without a named effect, the weighted genetic index
is used; with a named effect, selection is based on that effect.

Phenotypic selection uses previously stored phenotypic values and does not
re-evaluate the population. Without a named effect, the weighted phenotypic
index is used; with a named effect, its effect-specific phenotype is used.
The argument \code{x} counts distinct values, so ties are retained together.
Continuous phenotypic values will ordinarily produce fewer ties than many
genetic scoring schemes.
}

\seealso{
\code{phenotype.population}, \code{population.sort},
\code{select.n.best}
}
