\name{population.info.set}
\alias{population.info.set}
\title{Set stored information for an individual}

\description{
Sets the character information stored for an individual in a simulation population.
}

\usage{
population.info.set(PopName, ind, info)
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{ind}{Individual index or indices.}
   \item{info}{Information value.}
}
\value{
For a valid individual index, an invisible implementation-level \code{list} returned by the compiled setter. For an invalid index, \code{NULL}. The main effect is updating the stored individual information.
}
