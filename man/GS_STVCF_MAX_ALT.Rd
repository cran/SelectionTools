\name{.GS_STVCF_MAX_ALT}
\alias{.GS_STVCF_MAX_ALT}
\title{Maximum Number of ALT Alleles in STvcf Data}
\description{
Internal constant giving the maximum number of ALT alleles supported for one
marker in an STvcf object.
}
\usage{
.GS_STVCF_MAX_ALT
}
\format{
An integer scalar with value \code{253}.
}
\details{
Together with the REF allele, at most 254 biological allele states are
supported for one marker.  Missing data use a separate internal state.
This constant is intended for internal use and should not be modified.
}
\value{
An integer scalar with value \code{253}, giving the maximum number of ALT
alleles supported for one marker in an STvcf object. Together with the REF
allele this permits 254 biological allele states.
}
\keyword{internal}
