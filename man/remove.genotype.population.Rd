\name{remove.genotype.population}
\alias{remove.genotype.population}
\title{Removes cached genotype information from a population}

\description{
Clears cached marker-genotype data for one or more simulation populations.
}

\usage{
remove.genotype.population(PopNames)
}

\arguments{
\item{PopNames}{Population name or population specification.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of removing cached genotype information.
}

\details{
The raw chromosome-segment representation is retained. Existing genetic
values are retained. Stored phenotypic values are discarded, because
phenotypic values are tied to the current
genetic/genotype state. Genotype-dependent operations must regenerate marker
genotypes from the current linkage map and chromosome data.
}
