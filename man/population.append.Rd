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
\item{NameP2}{Character string naming the source population, which is not
modified.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of appending a copy of one simulation population to another.
}

\details{
Raw chromosome data are copied for every appended individual. If the
destination is initially empty, complete genotype, genetic-value, and
phenotypic-value categories are copied when present in the source. For a
nonempty destination, a derived category is retained only when it is complete
in both populations. Otherwise that category is cleared in the result.
}
