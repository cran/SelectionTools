\name{population.list}
\alias{population.list}
\title{List a population or population data}

\description{
Lists the currently supported simulation populations and their active sizes.
}

\usage{
population.list()
}
\value{
A \code{data.frame} with columns \code{PopName} and \code{count}, giving each current population name and its number of individuals.
}

\details{
The returned population registry is independent of SelectionTools marker-data-set names. The count associated with each name is the active number of individuals, not necessarily the allocated storage capacity.
}
