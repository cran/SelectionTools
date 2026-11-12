\name{population.copy}
\alias{population.copy}
\title{Copy a simulation population}

\description{
Copies all or a contiguous subset of a source population into a destination
population.
}

\usage{
population.copy(NameP1, NameP2, start=1, n=-1)
}

\arguments{
\item{NameP1}{Character string naming the destination population.}
\item{NameP2}{Character string naming the source population.}
\item{start}{One-based position of the first source individual to copy.}
\item{n}{Number of individuals to copy. A negative value requests the
remainder of the source population.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of copying a simulation population.
}

\details{
The source population is not modified. Raw chromosome data are copied.
Complete cached genotype, genetic-value, and phenotypic-value categories are
deep-copied when present in the source. Derived categories are maintained as
all-or-none population data; allocation failure cannot leave only a subset of
active destination individuals with a category.
}
