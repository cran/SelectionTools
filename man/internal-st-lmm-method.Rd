\name{.st_lmm_method}
\alias{.st_lmm_method}
\title{Dispatches compiled Mixed-Model Methods}
\description{
Dispatches mixed-model accessors, inference and display operations to the
compiled implementation.
}
\usage{
.st_lmm_method(object, action, ...)
}
\arguments{
\item{object}{A fitted \code{st.lmm} object, or a result value for an
operation that adds, removes or prints result display classes.}
\item{action}{A character string identifying the operation to perform.}
\item{...}{Named arguments required by the selected operation.}
}
\value{
The result of the selected operation, such as coefficients, a covariance
matrix, fitted values, a test table or a model summary. Display-class
operations return the corresponding result value. Printing operations
return the supplied object.
}
\keyword{internal}
