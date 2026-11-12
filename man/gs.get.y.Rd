\name{gs.get.y}
\alias{gs.get.y}
\title{Returns the genomic-selection response vector for a data set}

\description{
Returns the genomic-selection response vector for a data set.
}

\usage{
gs.get.y(data.set = "default", names = FALSE)
}

\arguments{
   \item{data.set}{Name of the SelectionTools data set.}
   \item{names}{Logical. Attach individual identifiers as vector names.}
}
\value{
Invisibly returns the current numeric response vector \code{y} in internal
individual order, including missing responses. With \code{names = TRUE},
each element is named with its individual identifier. The default returns
an unnamed vector. Returns \code{NULL} if the response cannot be obtained
or an argument is invalid.
}
