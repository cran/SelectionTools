\name{info}
\alias{info}
\title{Prints or handles an information message}

\description{
Prints or handles an information message.
}

\usage{
info(lev, msg)
}

\arguments{
   \item{lev}{Information or verbosity level.}
   \item{msg}{Message text.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of issuing the requested information message.
}
