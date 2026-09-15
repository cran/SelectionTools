\name{population.rename}
\alias{population.rename}
\title{Rename a simulation population}

\description{
Renames a registered simulation population.
}

\usage{
population.rename(OldName, NewName)
}

\arguments{
   \item{OldName}{Function argument used by this operation.}
   \item{NewName}{New population name. The name must contain at least two and
fewer than 256 characters.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of renaming a simulation population.
}
