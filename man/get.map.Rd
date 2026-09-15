\name{get.map}
\alias{get.map}
\title{get.map()}

\description{
Returns the currently loaded simulation linkage map.
}

\usage{
get.map()
}

\details{
The result describes the simulation map points in their active access order, including chromosome, map position, locus name, and class. Retrieving the map does not modify simulation populations or marker-data sets.
}

\value{
Returns the same \code{data.frame} as \code{linkage.map.get()}, with columns \code{chrom}, \code{pos}, \code{name}, and \code{class}.
}
