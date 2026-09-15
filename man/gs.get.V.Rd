\name{gs.get.V}
\alias{gs.get.V}
\title{Return the genomic-selection variance or covariance matrix for a data set}

\description{
Return the genomic-selection variance or covariance matrix for a data set.
}

\usage{
gs.get.V(data.set = "default")
}

\arguments{
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns the current numeric variance/covariance \code{matrix} \code{V}; returns \code{NULL} if it cannot be obtained from the current data set.
}
