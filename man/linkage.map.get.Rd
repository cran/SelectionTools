\name{linkage.map.get}
\alias{linkage.map.get}
\title{Get a linkage map used by the simulation routines}

\description{
Returns the currently active simulation linkage map.
}

\usage{
linkage.map.get()
}

\details{
The result contains the map points used by the simulation backend, including chromosome, position, locus name, and class, in the current simulation map order.
}

\value{
A \code{data.frame} with columns \code{chrom}, \code{pos}, \code{name}, and \code{class}, describing the active simulation linkage map.
}
