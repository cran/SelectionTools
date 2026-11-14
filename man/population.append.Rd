\name{population.append}
\alias{population.append}
\title{Appends a copy of a simulation population}

\description{
Appends a copy of a source population to the end of a destination population.
}

\usage{
population.append(NameP1, NameP2)
}

\arguments{
\item{NameP1}{Character string naming the recipient population.}
\item{NameP2}{Character string naming the source population. Its individuals are retained; mixed derived caches are cleared silently.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of appending a copy of one simulation population to another.
}

\details{
Raw chromosome data are copied for every appended individual. Cached marker
genotypes, genetic values and phenotypes are cleared throughout the resulting
destination population, including when it was initially empty. The source
population is retained. A mixed source or destination cache causes all three
caches in that population to be cleared silently.
}
