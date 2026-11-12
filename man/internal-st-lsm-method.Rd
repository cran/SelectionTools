\name{.st_lsm_method}
\alias{.st_lsm_method}
\title{Dispatches compiled least-squares mean Methods}
\description{
Dispatches Least-squares mean summaries, coefficient and covariance accessors,
hypothesis tests and display operations to the compiled implementation.
}
\usage{
.st_lsm_method(object, action, ...)
}
\arguments{
\item{object}{An \code{st.lsm} object, or a \code{summary.st.lsm} data frame
for printing a stored summary.}
\item{action}{A character string identifying the operation to perform.}
\item{...}{Named arguments required by the selected operation.}
}
\value{
The result of the selected operation, such as a summary or test table,
a coefficient vector, a covariance matrix or an updated Least-squares mean
object. Printing operations return the supplied object.
}
\keyword{internal}
