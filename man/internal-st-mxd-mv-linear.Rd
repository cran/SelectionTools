\name{.st_mxd_mv_linear}
\alias{.st_mxd_mv_linear}
\title{Evaluates Linear Functions of compiled Fixed Effects}
\description{
Calculates linear functions of fitted fixed effects, their uncertainty
and their estimability from a compiled mixed-model fit.
}
\usage{
.st_mxd_mv_linear(object, L, tol = 1e-7)
}
\arguments{
\item{object}{A compiled matrix-engine fit, or a fitted model containing that
fit in its \code{engine} component.}
\item{L}{A finite numeric matrix with one column per fixed coefficient and
one row per linear function. A vector specifies a single function.
Coefficient names, when supplied, are used to align the columns.}
\item{tol}{A positive numerical tolerance less than one for checking
estimability.}
}
\value{
A list with components \code{estimate}, \code{SE}, \code{vcov},
\code{estimable} and the aligned matrix \code{L}. Non-estimable functions
have missing estimates and standard errors; their covariance rows and
columns are also missing.
}
\keyword{internal}
