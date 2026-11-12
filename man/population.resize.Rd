\name{population.resize}
\alias{population.resize}
\title{Resize a population or population data}

\description{
Changes the allocated size of a simulation population.
}

\usage{
population.resize(PopName, newSize)
}

\arguments{
\item{PopName}{Character string naming the population to resize.}
   \item{newSize}{Requested new population allocation.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of resizing a simulation population.
}

\details{
A successful change of population allocation invalidates cached marker genotypes and evaluated genetic values because the active/storage relationship has changed. If no actual resize occurs, these derived categories are not cleared. If an attempted growth fails and the original allocation remains unchanged, the existing population and its derived data are preserved.
}
