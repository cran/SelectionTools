\name{population.plabsim.load}
\alias{population.plabsim.load}
\title{Loads a population from PLABSIM-style data}

\description{
Loads a simulation population from the PLABSIM-style file representation and restores stored individual information.
}

\usage{
population.plabsim.load(file, PopName=NA)
}

\arguments{
   \item{file}{File name or file connection.}
   \item{PopName}{Name of the population.}
}
\value{
\code{NULL}. The function is called for its side effect of loading a population from the PLABSIM-style representation and restoring stored individual information.
}
