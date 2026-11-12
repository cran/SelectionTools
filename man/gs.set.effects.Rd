\name{gs.set.effects}
\alias{gs.set.effects}
\title{Sets effects used by the genomic-selection routines}

\description{
Sets effects used by the genomic-selection routines.
}

\usage{
gs.set.effects(eff, filename ="seteffects", data.set ="default")
}

\arguments{
   \item{eff}{Function argument used by this operation.}
   \item{filename}{Character string retained for compatibility. The effects are transferred in memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of setting marker effects in the genomic-selection data set. Validation failures may return \code{NULL}.
}
