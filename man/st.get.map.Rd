\name{st.get.map}
\alias{st.get.map}
\title{Return the linkage map for a SelectionTools data set}

\description{
Return the linkage map for a SelectionTools data set.
}

\usage{
st.get.map(filename = "map", data.set = "default")
}

\arguments{
   \item{filename}{Character string retained for compatibility and used as part of the name of the internal session-temporary transfer file. No output file is retained.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{Pos}, \code{Name}, and \code{Class} describing the map of \code{data.set}; returns \code{NULL} if no map can be returned.
}
