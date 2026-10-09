\name{st.read.map.df}
\alias{st.read.map.df}
\title{Load a marker map from a data frame or matrix}
\description{
  Attaches a linkage map supplied in R to an existing marker data set,
  without writing temporary files.
}
\usage{
  st.read.map.df(map, m.stretch = 1, data.set = "default")
}
\arguments{
  \item{map}{A non-empty data frame or matrix with exactly three columns,
    interpreted in order as marker name, chromosome, and position (MCP).
    Column names are ignored. Marker names must be unique, nonempty, and contain
    no whitespace. Chromosomes must be positive integers and positions finite,
    non-negative numbers in the same units used by \code{st.read.map()}.}
  \item{m.stretch}{Positive finite multiplier applied to map positions.}
  \item{data.set}{Name of an existing marker data set.}
}
\details{
  Only markers shared by the map and marker data are retained. Genotypes are
  reordered by chromosome and position together with their marker names.
  At least one marker must match. Equal positions are separated by small
  internal offsets using the same rules as \code{st.read.map()}.
  
  Individual order, performance data, and the phase state are preserved.
  Marker statistics are recalculated; previous marker-dependent model results
  are cleared. Invalid input leaves the existing data set unchanged.
  This function attaches a marker-data map; it does not define the simulation genome.
}
\value{
  Returns \code{NULL} invisibly. The named marker data set is updated.
}
