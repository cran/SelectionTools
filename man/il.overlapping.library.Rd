\name{il.overlapping.library}
\alias{il.overlapping.library}
\title{Constructs or evaluates introgression-line libraries}

\description{
Constructs or evaluates introgression-line libraries.
}

\usage{
il.overlapping.library(n.c = 5, l.c = 100, s.l = 20, data.set = "default")
}

\arguments{
   \item{n.c}{Number or size used by the operation.}
   \item{l.c}{Function argument used by this operation.}
   \item{s.l}{Function argument used by this operation.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
An integer cleanup status returned by the final \code{unlink()} call after construction of the overlapping introgression-line library. The main result is the population/library created by the function.
}
