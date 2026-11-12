\name{gs.get.im}
\alias{gs.get.im}
\title{Returns the genomic-selection information matrix for a data set}

\description{
Returns the genomic-selection information matrix for a data set.
}

\usage{
gs.get.im(data.set = "default")
}

\arguments{
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
An invisible named integer vector with components \code{i}, the number of individuals, and \code{m}, the number of markers in the requested data set.
}
