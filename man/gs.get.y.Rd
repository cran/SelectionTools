\name{gs.get.y}
\alias{gs.get.y}
\title{Return the genomic-selection response vector for a data set}

\description{
Return the genomic-selection response vector for a data set.
}

\usage{
gs.get.y(data.set = "default")
}

\arguments{
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns the current numeric response vector \code{y}; returns \code{NULL} if it cannot be obtained from the current data set.
}
