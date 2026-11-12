\name{gs.get.Z}
\alias{gs.get.Z}
\title{Returns the genomic-selection marker or design matrix for a data set}

\description{
Returns the genomic-selection marker or design matrix for a data set.
}

\usage{
gs.get.Z(data.set = "default", names = FALSE)
}

\arguments{
   \item{data.set}{Name of the SelectionTools data set.}
   \item{names}{Logical. Attach individual identifiers as matrix row names.}
}
\value{
Invisibly returns the current numeric marker/design \code{matrix} \code{Z}
in internal individual order. With \code{names = TRUE}, row names are the
individual identifiers, in the same order as \code{gs.get.y(data.set, names = TRUE)}.
Column names are not added. The default returns an unnamed matrix.
Returns \code{NULL} if the matrix cannot be obtained or an argument is invalid.
}
\details{
\code{gs.make.A} carries the row names of \code{Z} to both dimensions of
the relationship matrix, so no separate name assignment is needed.
The option changes only the returned matrix attributes; marker values,
individual order and the internal data set are unchanged.
}
\seealso{\code{\link{gs.get.y}}, \code{\link{gs.make.A}}}
\examples{
data("v-tropmaize-vcf", package = "SelectionTools")
st.load.vcf.data(v.tropmaize.vcf, data.set = "maize.names")
Z <- gs.get.Z("maize.names", names = TRUE)
A <- gs.make.A(Z)
head(rownames(Z))
identical(rownames(A), rownames(Z))
}
