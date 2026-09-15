\name{gs.set.effects}
\alias{gs.set.effects}
\title{Set effects used by the genomic-selection routines}

\description{
Set effects used by the genomic-selection routines.
}

\usage{
gs.set.effects(eff, filename ="seteffects", data.set ="default")
}

\arguments{
   \item{eff}{Function argument used by this operation.}
   \item{filename}{Character string retained for compatibility and used as part of the name of the internal session-temporary file used to pass the effects to compiled code. The file is removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of setting marker effects in the genomic-selection data set. Validation failures may return \code{NULL}.
}
