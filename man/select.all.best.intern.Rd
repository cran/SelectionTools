% file select.all.best.intern.Rd
% copyright (C) 2001 Frisch/Maurer
%
\name{select.all.best.intern}
\alias{select.all.best.intern}
\title{Determines the best value classes in a sorted population}

\description{
Determines how many individuals belong to the requested number of best
value classes in an already sorted simulation population.
}

\usage{
select.all.best.intern(NamePop, n=1, selection.criterion="G", effect=NULL)
}

\arguments{
\item{NamePop}{Character string naming the already sorted population.}
\item{n}{Number of distinct value classes to retain.}
\item{selection.criterion}{Character string. \code{"G"}, the default,
uses genetic values; \code{"P"} uses phenotypic values.}
\item{effect}{Optional character string naming the effect whose genetic or
phenotypic value is used. If omitted, element 0 of the selected value array
is used.}
}

\details{
The function does not sort the population. It is used internally after
\code{population.sort}. The first \code{n} distinct values in the current
population order are treated as the best classes, and ties within a class
are retained together.

With \code{selection.criterion="G"}, genetic values are inspected. With
\code{selection.criterion="P"}, previously generated phenotypic values are
inspected. If \code{effect} is omitted, element 0 of the corresponding value
array is used; otherwise the named effect-specific value is used.
}

\value{
On success, a one-row \code{data.frame} with columns \code{n}, \code{minscore}, and \code{maxscore}, giving the number of individuals in the selected class(es) and the corresponding score limits. Returns \code{NULL} invisibly if the compiled selection step fails.
}

\seealso{
\code{select.all.best}, \code{population.sort}
}
