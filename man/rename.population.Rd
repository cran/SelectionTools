\name{rename.population}
\alias{rename.population}
\title{rename.population()}

\description{
Changes the registered name of a simulation population without changing its individuals.
}

\usage{
rename.population(OldName, NewName)
}

\arguments{
   \item{OldName}{Name of the population to be renamed}
   \item{NewName}{New name of the population}
}
\value{
Returns the same invisible implementation-level \code{list} as \code{population.rename()}. The function is a compatibility wrapper and the meaningful result is the delegated population operation.
}

\details{
The population object and its chromosome, genotype, evaluation, and information data remain the same; only the population identifier changes.
}

