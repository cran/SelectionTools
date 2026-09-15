\name{gd.splitdot}
\alias{gd.splitdot}
\title{Split or transform marker labels used by genetic-diversity routines}

\description{
Split or transform marker labels used by genetic-diversity routines.
}

\usage{
gd.splitdot(m.a)
}

\arguments{
   \item{m.a}{Marker or matrix input.}
}
\value{
A character \code{matrix} with two columns and one row for each input string. The first column contains the part before the first dot (or underscore if no dot is present), and the second column contains the remaining part; if no separator is present the second field is empty.
}
