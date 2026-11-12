\name{.st_mxd_mv}
\alias{.st_mxd_mv}
\title{Fits a Gaussian Mixed Model from Matrix Inputs}
\description{
Passes response, design and covariance-block specifications to the compiled
Gaussian mixed-model fitting engine.
}
\usage{
.st_mxd_mv(y, X, random = list(), residual = NULL,
           method = "REML", control = list())
}
\arguments{
\item{y}{A numeric response vector or matrix. Matrix columns represent
traits. Missing response cells are allowed; at least two observed cells
are required.}
\item{X}{A finite numeric fixed-effects design matrix with one row per
response cell, including missing cells. For a matrix response, rows
follow the column-wise ordering of \code{y}.}
\item{random}{A list of random-effect covariance-block specifications.}
\item{residual}{A list of residual covariance-block specifications, or
\code{NULL} for the default residual structure.}
\item{method}{A character string, either \code{"REML"} or \code{"ML"}.}
\item{control}{A named list of controls for the compiled fitting engine.}
}
\details{
Random and residual blocks describe covariance in the same response-cell
order as the design matrix. With \code{residual = NULL}, a vector response
receives independent residuals, while a matrix response receives an
unstructured covariance among traits within each response row.
}
\value{
A named list containing fixed-effect estimates and their covariance,
random-effect predictions, covariance-component estimates, fitted values,
residuals, likelihood information and convergence diagnostics.
}
\keyword{internal}
