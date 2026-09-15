% copyright (C) 2001 Frisch/Maurer
%
\name{gs.return.effects}
\alias{gs.return.effects}
\title{gs.return.effects()}
\description{
Returns the marker effects currently stored in a SelectionTools data set. The
underlying routine writes the effects to an intermediate file; the R wrapper
reads that file into R and deletes it before returning.
}
\usage{
gs.return.effects(out.filename ="effects", data.set ="default")
}
\arguments{
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{data.set}{character. Name of the data set that is used by the function.}
}

%\details{
%
%}


\value{
A \code{data.frame} containing the marker-effect table. If the underlying routine reports an error, \code{NULL} is returned invisibly. The intermediate file used to transfer the table is deleted after reading.
}

\note{
This function provides access to effects stored by genomic-effect estimation
functions such as \code{gs.esteff.rr()}, \code{gs.esteff.rmlc()}, and
\code{gs.esteff.rmlv()}. Their \code{auxfiles} argument defaults to
\code{FALSE}; \code{gs.return.effects()} can be used to retrieve the stored
effects independently of that output option.
}
%
