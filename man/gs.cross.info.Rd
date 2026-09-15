\name{gs.cross.info}
\alias{gs.cross.info}
\title{Return or write information for evaluated crosses}

\description{
Return or write information for evaluated crosses.
}

\usage{
gs.cross.info(bestn = 0, sortby = "index", out.filename ="cross.info",
              data.set ="default")
}

\arguments{
   \item{bestn}{Number of best entries to retain.}
   \item{sortby}{Criterion used for sorting.}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} containing the evaluated-cross information read from the generated result when the compiled operation succeeds; otherwise \code{NULL}.
}
