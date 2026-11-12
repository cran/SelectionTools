\name{.st_lmm_plain}
\alias{.st_lmm_plain}
\title{Restore Original Classes of Mixed-Model Results}
\description{
Removes mixed-model result display classes and restores the original
explicit classes recorded when the results were wrapped.
}
\usage{
.st_lmm_plain(x, recursive = FALSE)
}
\arguments{
\item{x}{An R object, possibly carrying the \code{st.lmm.result} class.}
\item{recursive}{Logical. If \code{TRUE}, also processes list components
recursively.}
}
\value{
The object with applicable result display classes and their bookkeeping
attributes removed. Original classes, values and other attributes are
retained. Objects requiring no conversion are returned unchanged.
}
\keyword{internal}
