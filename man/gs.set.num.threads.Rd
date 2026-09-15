\name{gs.set.num.threads}
\alias{gs.set.num.threads}
\title{Set num threads used by the genomic-selection routines}

\description{
Set num threads used by the genomic-selection routines.
}

\usage{
gs.set.num.threads(active)
}

\arguments{
   \item{active}{Function argument used by this operation.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of setting the number of genomic-selection computational threads.
}
