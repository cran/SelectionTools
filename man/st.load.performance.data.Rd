\name{st.load.performance.data}
\alias{st.load.performance.data}
\title{Load Performance Data from an R Data Set}
\description{
Matches performance records supplied in an R object to the individuals of an
existing SelectionTools marker data set without intermediate text files.
}
\usage{
st.load.performance.data(performance.data, data.set="default")
}
\arguments{
\item{performance.data}{A data frame or matrix containing individual
identifiers and numeric performance values.}
\item{data.set}{Character string naming the SelectionTools marker data set.}
}
\details{
Marker data must already be present in \code{data.set}. Three input-column
conventions are accepted. Columns named \code{individual} and
\code{performance} are used when both are present; otherwise columns
\code{ii} and \code{yy} are used when both are present; otherwise an object
with exactly two columns is interpreted as identifier followed by performance.

Performance individual IDs must contain only ASCII letters and digits and
fewer than 256 bytes. The same name convention is used for STvcf individual
identifiers. Invalid IDs are rejected rather than being silently ignored as
unmatched records.

Records are matched by individual identifier. Performance records whose valid
IDs do not occur in the current marker data are ignored. Marker-data
individuals with no matched performance record are removed. If an identifier
occurs more than once in the performance input, the last matched value is
used, following the existing performance-input convention. At least one
record must match.

Performance values must be finite numeric values after conversion. The
resulting marker data set retains the current marker names and linkage map,
subsets both allele arrays to individuals with performance, recalculates
marker statistics for the retained individuals, and stores the matched
performance vector. The operation is transactional.
}
\value{
\code{NULL}, returned invisibly. The function is called for its side effect of loading performance data into the named SelectionTools data set.
}
