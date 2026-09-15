\name{gs.get.Z}
\alias{gs.get.Z}
\title{Return the genomic-selection marker or design matrix for a data set}

\description{
Return the genomic-selection marker or design matrix for a data set.
}

\usage{
gs.get.Z(data.set = "default")
}

\arguments{
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns the current numeric marker/design \code{matrix} \code{Z}; returns \code{NULL} if it cannot be obtained from the current data set.
}
