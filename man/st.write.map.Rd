\name{st.write.map}
\alias{st.write.map}
\title{Writes linkage-map information for a SelectionTools data set}

\description{
Writes linkage-map information for a SelectionTools data set.
}

\usage{
st.write.map(mfilename, auxfiles = TRUE, data.set = "default")
}

\arguments{
   \item{mfilename}{Required base name for the map file. The file is written with extension \file{.mmp}.}
   \item{auxfiles}{Logical. If \code{TRUE}, map information is also returned as a data frame assembled in memory. The map file itself is written regardless of this argument.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{Pos}, \code{Name}, and \code{Class} when \code{auxfiles = TRUE} and map output succeeds; otherwise \code{NULL}. The map file is also written as requested.
}
