\name{set.population.info}
\alias{set.population.info}
\title{Set population information used by the simulation routines}

\description{
Set population information used by the simulation routines.
}

\usage{
set.population.info(PopName, ind, info)
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{ind}{Individual index or indices.}
   \item{info}{Information value.}
}
\value{
Returns the same invisible implementation-level \code{list} or \code{NULL} as \code{population.info.set()}, depending on whether the individual index is valid.
}
