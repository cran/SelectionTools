\name{list.populations}
\alias{list.populations}
\title{list.populations}

\description{
Returns the names and active sizes of the currently registered simulation populations.
}

\usage{
list.populations()
}
\value{
Returns the same \code{data.frame} as \code{population.list()}, with columns \code{PopName} and \code{count}.
}

\details{
The result reflects the simulation population registry and is independent of SelectionTools marker-data-set names.
}

