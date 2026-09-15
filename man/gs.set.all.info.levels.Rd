\name{gs.set.all.info.levels}
\alias{gs.set.all.info.levels}
\title{Set all info levels used by the genomic-selection routines}

\description{
Set all info levels used by the genomic-selection routines.
}

\usage{
gs.set.all.info.levels(level)
}

\arguments{
   \item{level}{Function argument used by this operation.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of setting all genomic-selection information levels.
}
