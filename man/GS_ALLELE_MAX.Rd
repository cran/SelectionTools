\name{.GS_ALLELE_MAX}
\alias{.GS_ALLELE_MAX}
\title{Maximum Allele Code Used by SelectionTools}
\description{
Internal constant giving the largest biological allele code supported by
SelectionTools.
}
\usage{
.GS_ALLELE_MAX
}
\format{
An integer scalar with value \code{254}.
}
\details{
The value is used when validating allele codes in the R interface.  The
separate internal missing-value code is not a biological allele code.
This constant is intended for internal use and should not be modified.
}
\value{
An integer scalar with value \code{254}, giving the largest biological allele
code supported by the SelectionTools R interface. The separate missing-value
code is not included in this maximum.
}
\keyword{internal}
