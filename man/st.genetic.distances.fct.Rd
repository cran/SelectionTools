\name{st.genetic.distances.fct}
\alias{st.genetic.distances.fct}
\title{Calculate genetic distances from SelectionTools marker data}

\description{
Calculate genetic distances from SelectionTools marker data.
}

\usage{
st.genetic.distances.fct(measure = "mrd", format = "l", split = 1,
                         filename = NULL, auxfiles = FALSE,
                         data.set = "default")
}

\arguments{
   \item{measure}{Distance or similarity measure.}
   \item{format}{Input or output format.}
   \item{split}{Split or partition control.}
   \item{filename}{Character string giving the base name of the genetic-distance file. Required when \code{auxfiles = TRUE}; the file is written with extension \file{.gdi} in the output directory.}
   \item{auxfiles}{Logical. If \code{TRUE}, a persistent output file with the genetic distances is generated and \code{filename} must be supplied.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns either an object of class \code{"dist"} when \code{format = "m"}, representing pairwise genetic distances, or a long-format \code{data.frame} with columns \code{OTU1}, \code{OTU2}, and \code{Measure} when \code{format = "l"}.
}
