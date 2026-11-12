\name{gs.cross.info}
\alias{gs.cross.info}
\title{Information for Evaluated Crosses}

\description{
Returns cross summaries with full numerical precision.
}

\usage{
gs.cross.info(bestn = 0, sortby = "index", out.filename = "cross.info",
              data.set = "default")
}

\arguments{
   \item{bestn}{Nonnegative whole number of entries to retain after sorting.
   A value of 0 returns all crosses. Values larger than the number of crosses
   return all available entries.}
   \item{sortby}{Character string naming the sorting criterion:
   \code{"index"}, \code{"gd"}, \code{"mu"}, \code{"mi"}, \code{"ma"},
   \code{"va"} or \code{"es"}. The criterion must already have been
   evaluated unless \code{"index"} is used.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} containing the cross summaries.
Columns \code{P1Name} and \code{P2Name} are character strings, so
numeric-looking individual identifiers are retained unchanged.
Columns \code{P1No} and \code{P2No} contain parent indices. The numeric
columns \code{gd}, \code{mu}, \code{mi}, \code{ma}, \code{va} and \code{es}
contain genetic distance, expected mean, lower and upper marker-effect limits,
segregation variance and expected superior progeny value, respectively.
Unevaluated criteria are \code{NA}; an evaluated zero remains zero.
The segregation standard deviation is \code{sqrt(va)}.
}
\details{
The \code{"index"} order follows the parent indices. The \code{"mi"}
criterion is sorted in increasing order; the other numerical criteria are
sorted in decreasing order. Ties follow the parent indices.

Numeric results are returned directly from the stored calculations without
rounding through formatted text. Display precision is controlled by R's
printing options and does not alter the returned values.

Changing the marker data, map or marker effects invalidates affected cross
summaries. Subsequent evaluations use the current data and effects.
Invalid arguments and unavailable sorting criteria produce an error.
}
\seealso{
\code{gs.cross.eval.mu}, \code{gs.cross.eval.va}, \code{gs.cross.eval.es}
}
