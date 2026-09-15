\name{load.linkage.map}
\alias{load.linkage.map}
\title{load.linkage.map()}

\description{
Loads the simulation linkage map from a text file.
}

\usage{
load.linkage.map(file, disperse=TRUE, file.disperse=NA, disperse.factor=100)
}

\arguments{
\item{file}{Character string naming the simulation linkage-map input file.}
   \item{disperse}{Logical. If true, coincident locus positions are deterministically separated during map preparation.}
   \item{file.disperse}{Optional character string naming a file to which the dispersed map is written. An omitted value does not create an additional file.}
   \item{disperse.factor}{Finite numeric factor controlling the spacing used when coincident positions are dispersed.}
}
\value{
Returns the same invisible status \code{list} as \code{linkage.map.load()}, including its \code{retval} component.
}

\details{
The simulation linkage map is distinct from the marker-data map managed by \code{st.read.map}. It defines the loci available to simulation-population genotyping, recombination, and evaluation.

Map records contain chromosome, position, locus name, and class. The loader prepares the loci in chromosome/position order and can disperse coincident map points according to its arguments. Only one simulation linkage map is active at a time.

The function \code{linkage.map.load} provides the same linkage-map loading operation.
}
