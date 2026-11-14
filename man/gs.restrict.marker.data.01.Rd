\name{gs.restrict.marker.data.01}
\alias{gs.restrict.marker.data.01}
\title{Restricts marker data using the package genomic-selection interface}

\description{
Restricts marker data using the package genomic-selection interface.
}

\usage{
gs.restrict.marker.data.01(outfile, traitfile = NULL, NoAll.MAX = 999,
                           MaMis.MAX = 1, ExHet.MIN = 0, InMis.MAX = 1,
                           data.set = "default")
}

\arguments{
   \item{outfile}{Function argument used by this operation.}
   \item{traitfile}{Function argument used by this operation.}
   \item{NoAll.MAX}{Function argument used by this operation.}
   \item{MaMis.MAX}{Function argument used by this operation.}
   \item{ExHet.MIN}{Function argument used by this operation.}
   \item{InMis.MAX}{Function argument used by this operation.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of restricting the genomic-selection marker data. Validation failures may return \code{NULL}.
}
