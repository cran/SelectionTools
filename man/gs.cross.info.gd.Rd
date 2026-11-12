\name{gs.cross.info.gd}
\alias{gs.cross.info.gd}
\title{Returns or write genetic-distance information for crosses}

\description{
Returns or write genetic-distance information for crosses.
}

\usage{
gs.cross.info.gd(out.filename ="cross.info.gd", data.set ="default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} containing the genetic-distance cross information read from the generated result when the compiled operation succeeds; otherwise \code{NULL}.
}
