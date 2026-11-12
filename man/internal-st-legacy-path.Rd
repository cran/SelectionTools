\name{.st.legacy.path}
\alias{.st.legacy.path}
\title{Resolve a Legacy Marker-Data Path}
\description{
Resolves a file name using the configured input or output directory.
}
\usage{
.st.legacy.path(filename, output = FALSE)
}
\arguments{
\item{filename}{The file name to resolve.}
\item{output}{A logical flag. \code{TRUE} selects the output directory
and \code{FALSE} selects the input directory.}
}
\value{
The resolved file path as a character string.
}
\keyword{internal}
