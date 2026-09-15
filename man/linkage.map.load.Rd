\name{linkage.map.load}
\alias{linkage.map.load}
\title{Load a linkage map used by the simulation routines}

\description{
Loads and prepares the linkage map used by the simulation backend.
}

\usage{
linkage.map.load(file, disperse=TRUE, file.disperse=NA, disperse.factor=100)
}

\arguments{
\item{file}{Character string naming the simulation map file.}
   \item{disperse}{Logical flag controlling dispersion of coincident map positions.}
   \item{file.disperse}{Optional file-related dispersion specification used by the existing map loader.}
   \item{disperse.factor}{Numeric factor controlling the dispersion scale used by the existing loader.}
}
\value{
An invisible \code{list} with a \code{retval} status component. Invalid input returns \code{list(retval = -2L)} invisibly; otherwise the list is returned by the compiled map loader. The main effect is loading the linkage map.
}

\details{
The simulation map supplies chromosome, position, locus-name, and class information used by population genotyping and recombination. The loader prepares map points in a valid access order and can disperse coincident positions when requested. The active simulation map is a simulation-wide object and is separate from a SelectionTools marker data set.
}
