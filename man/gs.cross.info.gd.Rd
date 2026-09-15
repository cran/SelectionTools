\name{gs.cross.info.gd}
\alias{gs.cross.info.gd}
\title{Return or write genetic-distance information for crosses}

\description{
Return or write genetic-distance information for crosses.
}

\usage{
gs.cross.info.gd(out.filename ="cross.info.gd", data.set ="default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} containing the genetic-distance cross information read from the generated result when the compiled operation succeeds; otherwise \code{NULL}.
}
