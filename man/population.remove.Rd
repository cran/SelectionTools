\name{population.remove}
\alias{population.remove}
\title{Removes a population or population data}

\description{
Removes one or more simulation populations.
}

\usage{
population.remove(PopNames)
}

\arguments{
   \item{PopNames}{Population name or population specification.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of removing simulation populations.
}

\details{
The named populations are deleted from the simulation population registry and
their owned chromosome, information, genotype, and evaluation storage is
returned to the package's compiled memory pool. Other simulation-wide definitions
such as the genome, linkage map, and effects are separate and are not removed
by this operation.

The pool retains backing memory for subsequent simulation populations. Its
default backing-block increment is 100 MB (100000000 bytes); a larger
contiguous request can require a multiple of that increment. Population
deletion therefore need not reduce the process's resident memory. The pool's
backing memory is released when the package's compiled library is unloaded.
}
