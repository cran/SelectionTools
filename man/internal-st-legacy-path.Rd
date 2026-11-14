\name{.st.legacy.path}
\alias{.st.legacy.path}
\alias{.st.restore.graphics}
\title{Path Resolution and Graphics Restoration}
\description{
Resolves a file name using the configured input or output directory and
restores saved graphics parameters after plotting.
}
\usage{
.st.legacy.path(filename, output = FALSE)
.st.restore.graphics(settings)
}
\arguments{
\item{filename}{The file name to resolve.}
\item{output}{A logical flag. \code{TRUE} selects the output directory
and \code{FALSE} selects the input directory.}
\item{settings}{The complete named list returned by
\code{graphics::par(no.readonly=TRUE)} for the active graphics device.}
}
\details{
The graphics helper restores saved layout, margins and text scaling without
directly setting the derived \code{pin} and \code{plt} values. This permits
restoration of the original margins on small devices. The snapshot belongs
to the active device; the helper does not open or close graphics devices.
}
\value{
\code{.st.legacy.path} returns the resolved file path as a character string.
\code{.st.restore.graphics} returns \code{NULL}, invisibly, after restoring
the graphics settings.
}
\keyword{internal}
