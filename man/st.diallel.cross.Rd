\name{st.diallel.cross}
\alias{st.diallel.cross}
\title{Identify Unordered Diallel Crosses}
\description{Constructs a common cross factor for the two parental roles in
a half diallel.}
\usage{
st.diallel.cross(mother, father)
}
\arguments{
\item{mother,father}{Equal-length vectors of parent identifiers. Factors
are converted to their labels. Missing either parent produces a missing
cross.}
}
\details{
Swapping the parental roles leaves the cross identifier unchanged.
Length-prefixed labels avoid collisions when parent names contain the
separator. The labels are identifiers, not pedigree expressions.

\code{random=~overlay(mother, father)} in \code{\link{st.lmm}} specifies
general combining ability. The cross factor enters as a separate random
term for specific combining ability.
The overlay uses one shared parental effect vector; each parent has loading
one and a self cross has loading two. Reciprocal and maternal effects must
be specified separately when required.
}
\value{A factor of canonical unordered parent pairs, with missing values
where either parent identifier is missing.}
\seealso{\code{\link{st.lmm}}, \code{\link{st.predict}}}
