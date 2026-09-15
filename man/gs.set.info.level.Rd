\name{gs.set.info.level}
\alias{gs.set.info.level}
\title{Set info level used by the genomic-selection routines}

\description{
Set info level used by the genomic-selection routines.
}

\usage{
gs.set.info.level(level, data.set ="default")
}

\arguments{
   \item{level}{Function argument used by this operation.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of setting the genomic-selection information level for the data set.
}
