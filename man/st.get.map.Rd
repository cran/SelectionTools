\name{st.get.map}
\alias{st.get.map}
\title{Returns the linkage map for a SelectionTools data set}

\description{
Returns the genetic map, or the stored physical map when no genetic map
is present, for a SelectionTools data set.
}

\usage{
st.get.map(filename = "map", data.set = "default")
}

\arguments{
   \item{filename}{Character string retained for compatibility. The map is returned directly from memory; this argument does not create a file.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\details{
A physical map imported from STvcf is stored independently of genetic
positions. When it is returned, \code{Class = "bp"} identifies the
base-pair units of \code{Pos}. Genetic-map returns retain their previous
interpretation. Only markers currently present in the dataset are returned.
}
\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{Pos}, \code{Name}, and \code{Class} describing the map of \code{data.set}; returns \code{NULL} if no map can be returned.
}
