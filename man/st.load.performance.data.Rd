\name{st.load.performance.data}
\alias{st.load.performance.data}
\alias{st.generate.phenotypes}
\title{Loads Performance Data from an R Data Set}
\description{
Matches performance records supplied in an R object to the individuals of an
existing SelectionTools marker data set without intermediate text files.
}
\usage{
st.load.performance.data(performance.data, data.set="default")

st.generate.phenotypes(data.set, values, h2 = 0.8, trait = "y")
}
\arguments{
\item{performance.data}{A data frame or matrix containing individual
identifiers and numeric performance values.}


\item{data.set}{Name of the destination marker or supported hybrid dataset.}
\item{values}{Finite numeric vector in dataset order, or a prediction data frame returned by \code{st.predict}. For a data frame, identifiers must match the target rows exactly.}
\item{h2}{Heritability, greater than zero and at most one.}
\item{trait}{Phenotype column for a hybrid dataset. For an ordinary marker dataset the stored response is \code{y}. A parental identifier column cannot be replaced.}

}
\details{
Marker data must already be present in \code{data.set}. Three input-column
conventions are accepted. Columns named \code{individual} and
\code{performance} are used when both are present; otherwise columns
\code{ii} and \code{yy} are used when both are present; otherwise an object
with exactly two columns is interpreted as identifier followed by performance.

Performance individual IDs must be non-empty, subject to the R string-size limit. Whitespace and control characters are rejected; ordinary punctuation
such as dots, underscores, and hyphens is permitted. The same name convention
is used for STvcf individual identifiers. Invalid IDs are rejected rather than
being silently ignored as unmatched records.

Records are matched by individual identifier. Performance records whose valid
IDs do not occur in the current marker data are ignored. Marker-data
individuals with no matched performance record are removed. If an identifier
occurs more than once in the performance input, the last matched value is
used, following the existing performance-input convention. At least one
record must match.

Performance values may be finite numeric values or \code{NA}. An explicit
\code{NA} record retains that individual and its marker genotypes, with a
missing value in the stored performance vector. This differs from an
individual without any matching performance row, which is removed.
The character token \code{"NA"} is also accepted; other nonnumeric values,
\code{NaN} and infinities are rejected. The
resulting marker data set retains the current marker names and linkage map,
subsets both allele arrays to individuals with a matching record, recalculates
marker statistics for the retained individuals, and stores the matched
performance vector. The operation is transactional.
}
\value{
\code{NULL}, returned invisibly. The function is called for its side effect of loading performance data into the named SelectionTools data set.
}

\examples{
data("v-tropmaize-vcf", package = "SelectionTools")
marker.frame <- st.STvcf.to.dataframe(v.tropmaize.vcf)[1:12, c(1:6, 7:14)]
vcf <- st.dataframe.to.STvcf(marker.frame)
st.load.vcf.data(vcf, data.set = "maize.example")
data("v-tropmaize-phe", package = "SelectionTools")
performance <- v.tropmaize.phe[match(vcf$individual,
                                   as.character(v.tropmaize.phe$ii)), ]
names(performance) <- c("individual", "performance")
st.load.performance.data(performance, data.set = "maize.example")
}

\seealso{
\code{\link{st.read.performance.data}} for reading performance
data from a local file.
}

\section{Phenotype generation from genetic predictions}{
Adds independent normal residuals to supplied genetic predictions and stores the resulting phenotype in a marker or supported hybrid dataset.

With sample variance \eqn{s_g^2} of the supplied predictions, residual variance is \eqn{s_e^2=s_g^2(1-h^2)/h^2}. Finite-sample realized heritability need not equal the specified value. At least two records are required; predictions must vary when \code{h2 < 1}. Uses the R random-number generator. Dataset identifiers and pedigree are preserved. Replacing the response invalidates affected active fitted state. This operation does not send values to simulation populations.

The generated numeric phenotype vector, invisibly. The destination response is updated.

}
