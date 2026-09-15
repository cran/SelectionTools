\name{define.genome}
\alias{define.genome}
\title{Define a genome from a linkage map and optional chromosome lengths}

\description{
Define a genome from a linkage map and optional chromosome lengths.
}

\usage{
define.genome(map, length.vec=NULL)
}

\arguments{
   \item{map}{Linkage map or map object.}
   \item{length.vec}{Vector of chromosome lengths. If a chromosome length is zero, it is replaced by 0.01 and a warning is issued.}
}
\value{
\code{NULL} on successful completion and on the documented validation exits. The function is called for its side effects of setting genome parameters and defining the linkage map.
}
