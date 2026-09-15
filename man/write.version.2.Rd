\name{write.version.2}
\alias{write.version.2}
\title{Write version information in the package-specific format}

\description{
Write version information in the package-specific format.
}

\usage{
write.version.2()
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of requesting/printing package version information.
}
