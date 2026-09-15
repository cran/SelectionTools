\name{linkage.map.save}
\alias{linkage.map.save}
\title{Save a linkage map used by the simulation routines}

\description{
Writes a simulation linkage map to a text file.
}

\usage{
linkage.map.save(file, dta.map=linkage.map.get())
}

\arguments{
\item{file}{Character string naming the output map file.}
   \item{dta.map}{Map data to write; by default the currently active simulation linkage map.}
}
\value{
\code{NULL}, returned invisibly by \code{write.table()}. The function is called for its side effect of writing the linkage map.
}

\details{
By default the currently active simulation linkage map is obtained and written. A supplied \code{dta.map} can instead be written. Saving does not alter the active map or any simulation population.
}
