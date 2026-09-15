\name{return.population}
\alias{return.population}
\title{return.population()}

\description{
Converts one or more simulation populations to the marker-incidence data-frame representation used by the genetic-distance routines.
}

\usage{
return.population(PopNames, missing.rm=TRUE)
}

\arguments{
\item{PopNames}{Character string containing one or more simulation population names separated by blanks.}
   \item{missing.rm}{Logical flag controlling removal or retention of missing marker information in the returned representation.}
}

\details{
The requested simulation populations are evaluated at the current marker loci and returned in an NTSys-style matrix representation suitable for the corresponding SelectionTools data utilities. Multiple population names can be supplied in one space-separated string.

The returned object is a representation of population marker data; it does not remove or modify the source populations. The \code{missing.rm} argument controls handling of missing marker information in the conversion.
}
\value{
A \code{data.frame} in the marker-incidence representation used by the genetic-distance routines. Returns \code{NULL} invisibly if the compiled population conversion cannot produce a result.
}
