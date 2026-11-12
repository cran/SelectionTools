\name{population.sample}
\alias{population.sample}
\title{Sample a simulation population}

\description{
Creates a destination population by random sampling from a source population.
}

\usage{
population.sample(NameP1, NameP2, size=-1, replace=FALSE)
}

\arguments{
\item{NameP1}{Character string naming the sampled destination population.}
\item{NameP2}{Character string naming the source population.}
\item{size}{Requested sample size. A negative value requests the source
population size.}
\item{replace}{Logical. If true, sampling is with replacement.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of sampling a simulation population.
}

\details{
The source population is not changed. Raw chromosome data are copied for the
sampled individuals. Complete cached genotype, genetic-value, and
phenotypic-value categories are deep-copied when present in the source.
Derived categories remain all-or-none for the destination population.
}
