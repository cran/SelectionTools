\name{gs.info}
\alias{gs.info}
\title{Prints or handles a genomic-selection information message}

\description{
Prints or handles a genomic-selection information message.
}

\usage{
gs.info(lev, msg)
}

\arguments{
   \item{lev}{Information or verbosity level.}
   \item{msg}{Message text.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of issuing the requested genomic-selection information message.
}
