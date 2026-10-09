% copyright (C) 2001 Frisch/Maurer
%
\name{st.data.set.remove}
\alias{st.data.set.remove}
\title{Remove a SelectionTools Marker Data Set}
\description{
Releases all native storage belonging to one named SelectionTools marker data set.
}
\usage{
st.data.set.remove(data.set)
}
\arguments{
\item{data.set}{One non-empty character string naming the data set to remove.}
}
\details{
For a user-defined data set, all marker, map, phenotype, model, haplotype, and other derived native storage is released and the data-set name is removed from the registry. For \code{"default"}, the contents are reset but the built-in name remains registered.

The function does not create a missing data set. If \code{data.set} is not registered, no state is changed.
}
\value{
Invisibly returns \code{TRUE} when the named data set was found and cleared or removed, and \code{FALSE} when it was not registered.
}
