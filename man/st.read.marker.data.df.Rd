\name{st.read.marker.data.df}
\alias{st.read.marker.data.df}
\title{Read marker data from a data frame into SelectionTools}

\description{
Read marker data from a data frame into SelectionTools.
}

\usage{
st.read.marker.data.df(markerdata, format = NULL, data.set = "default")
}

\arguments{
   \item{markerdata}{Function argument used by this operation.}
   \item{format}{Input or output format.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
\code{NULL}, returned invisibly. The function is called for its side effect of loading the supplied marker-data frame into the named SelectionTools data set.
}

\details{
Biological allele codes must be integers from 1 through 254; they need not be
consecutive within a marker. R missing allele values and the character tokens
\code{-1}, \code{--}, \code{-}, \code{.}, and \code{NA} are treated as missing.
The \code{-1} spelling is retained for compatibility with historical
SelectionTools marker data. If
either allele of a diploid genotype is missing, the complete genotype is
stored as missing. Repeated individual/marker records are aggregated
deterministically using the same rules as list-format marker input.
}
