% copyright (C) 2001 Frisch/Maurer
%
\name{st.copy.marker.data}
\alias{st.copy.marker.data}
\title{st.copy.marker.data()}
\description{
Creates an independent copy of the primary marker data, marker statistics, linkage-map representation, and phenotype vector of another SelectionTools data set.
}
\usage{
st.copy.marker.data(target.data.set, source.data.set)
}
\arguments{
\item{target.data.set}{Character string naming the data set that will receive the copy. Existing contents of this target are replaced.}
   \item{source.data.set}{Character string naming the marker data set to copy.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of copying marker data from the source data set to the target data set.
}

\details{
The target receives its own marker-name and individual-name arrays, allele matrices, recomputed marker statistics, map storage, map access pointers, marker-to-map indices, and phenotype vector when present. Map storage order and map access order are reconstructed so that the copy does not depend on pointers into the source data set.

The copy is intended as a copy of the marker-data state, not of every fitted model or temporary calculation attached to the source. Model matrices, fitted effects, p-values, cross-validation state, and other derived analysis objects are not implied by this operation.

The source data set is not modified. Existing contents of the target data set are replaced.
}


%
%
