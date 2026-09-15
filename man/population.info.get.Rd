\name{population.info.get}
\alias{population.info.get}
\title{Return stored information for an individual}

\description{
Returns the character information stored for an individual in a simulation population.
}

\usage{
population.info.get(PopName, ind)
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{ind}{Individual index or indices.}
}
\value{
A character value containing the stored information for the requested individual when the individual index is valid; otherwise \code{NULL}.
}
