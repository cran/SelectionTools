% copyright (C) 2001 Frisch/Maurer
%
\name{st.restrict.marker.data}
\alias{st.restrict.marker.data}
\title{st.restrict.marker.data()}
\description{
Restricts an existing marker data set by individual names, marker names, missingness, allele number, and expected heterozygosity.
}
\usage{
st.restrict.marker.data(ind.list = NULL, ind.file = NULL, mar.list = NULL,
                        mar.file = NULL, NoAll.MAX = 999, MaMis.MAX = 1,
                        ExHet.MIN = 0, InMis.MAX = 1, data.set = "default")
}
\arguments{
\item{ind.list}{Optional character vector of individual names to retain. Unknown names are ignored. Do not supply together with \code{ind.file}.}
   \item{ind.file}{Optional file containing individual names to retain. Relative names are read from the input directory. Do not supply together with \code{ind.list}.}
   \item{mar.list}{Optional character vector of marker names to retain. Unknown names are ignored. Do not supply together with \code{mar.file}.}
   \item{mar.file}{Optional file containing marker names to retain. Do not supply together with \code{mar.list}.}
   \item{NoAll.MAX}{Non-negative maximum number of non-missing alleles allowed per marker.}
   \item{MaMis.MAX}{Finite maximum marker missingness.}
   \item{ExHet.MIN}{Finite minimum expected heterozygosity.}
   \item{InMis.MAX}{Finite maximum individual missingness.}
   \item{data.set}{Character string naming the marker data set to restrict.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of restricting the marker data in the selected data set. Validation failures may return \code{NULL}.
}

\details{
The function first constructs individual and marker inclusion masks from the current marker statistics. Individual criteria and marker criteria are evaluated against the statistics that exist before the restriction starts. The filtering is therefore a single-pass operation, not an iterative filter-and-recalculate procedure. Calling the function repeatedly can be used when sequential recalculation between criteria is desired.

If an individual or marker selection list is supplied, names not found in the data set are silently ignored. Duplicate names in a selection list have no additional effect. The order of a selection list does not reorder the marker data: surviving individuals and markers remain in their original data-set order.

Marker filters use the number of non-missing alleles, marker missingness, and expected heterozygosity. Individual filtering uses individual missingness. All active criteria are combined. If every marker and every individual satisfies the criteria, the function is a true no-op and leaves all existing derived structures untouched.

A genuine restriction is transactional. A complete temporary marker matrix, phenotype vector, restricted map, and marker statistics are built before the current data set is replaced. If no individual or no marker would remain, or if construction fails, the existing data set is left unchanged.

If a map is present, it is restricted through the map access order and marker-to-data mapping, rather than by assuming that physical map storage has the same order as the marker matrix. This preserves valid map structures even when storage order and marker order differ. If markers are removed, marker-dependent model state is invalidated. If only individuals are removed, a complete and internally consistent marker-effect definition can be retained; incomplete or inconsistent effect state is not retained.
}
\note{
Restrictions are based on the marker statistics present at the start of each call. Repeating the function after statistics have been recalculated can therefore give a different result from applying several criteria simultaneously.
}
%
