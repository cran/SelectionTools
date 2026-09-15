\name{il.ideal.library}
\alias{il.ideal.library}
\title{Construct or evaluate introgression-line libraries}

\description{
Construct or evaluate introgression-line libraries.
}

\usage{
il.ideal.library(n.c = 5, l.c = 100, s.l = 20, data.set = "default")
}

\arguments{
   \item{n.c}{Number or size used by the operation.}
   \item{l.c}{Function argument used by this operation.}
   \item{s.l}{Function argument used by this operation.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
An integer cleanup status returned by the final \code{unlink()} call after construction of the ideal introgression-line library. The main result is the population/library created by the function.
}
