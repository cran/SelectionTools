\name{.st_lmm_result}
\alias{.st_lmm_result}
\title{Adds Display Classes to Mixed-Model Results}
\description{
Adds the mixed-model result display class to supported values, retaining
the original classes for subsequent conversion and computation.
}
\usage{
.st_lmm_result(x)
}
\arguments{
\item{x}{An R result object to prepare for mixed-model result printing.}
}
\details{
Supported values include numeric, integer, logical and character objects,
lists and formulas. Eligible list components are processed recursively;
data-frame columns are not wrapped separately. Already wrapped results
and fitted model, Least-squares mean and prediction objects are retained.
}
\value{
The result with \code{st.lmm.result} display classes added where applicable.
Values and storage types are retained, and original classes are recorded.
Unsupported objects and \code{NULL} are returned unchanged.
}
\keyword{internal}
