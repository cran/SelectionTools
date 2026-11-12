\name{.st.valid.identifier}
\alias{.st.valid.identifier}
\title{Validate SelectionTools Identifiers}
\description{
Checks character identifiers used by marker-data input functions.
}
\usage{
.st.valid.identifier(x, max.bytes = 255L)
}
\arguments{
\item{max.bytes}{Maximum permitted identifier size in bytes; \code{Inf} permits individual identifiers subject to the R string-size limit.}
\item{x}{A character vector of identifiers.}
}
\details{
An identifier must be non-missing and non-empty, and its byte length must
not exceed \code{max.bytes}. Its UTF-8 representation must contain neither bytes with
values at most 32 nor the byte with value 127. This excludes ASCII spaces
and control characters.
}
\value{
A logical vector with one element for each input identifier. A value of
\code{TRUE} indicates that the identifier satisfies these conditions.
}
\keyword{internal}
