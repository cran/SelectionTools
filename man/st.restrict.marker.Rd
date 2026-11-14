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
\item{ind.list}{Optional character vector of individual names to retain. Unknown names are ignored. This is mutually exclusive with \code{ind.file}.}
   \item{ind.file}{Optional file containing individual names to retain. Relative names are read from the input directory. This is mutually exclusive with \code{ind.list}.}
   \item{mar.list}{Optional character vector of marker names to retain. Unknown names are ignored. This is mutually exclusive with \code{mar.file}.}
   \item{mar.file}{Optional file containing marker names to retain. This is mutually exclusive with \code{mar.list}.}
   \item{NoAll.MAX}{Maximum number of alleles allowed per marker.}
   \item{MaMis.MAX}{Maximum proportion of missing genotypes per marker.}
   \item{ExHet.MIN}{Minimum expected heterozygosity.}
   \item{InMis.MAX}{Maximum proportion of missing genotypes per individual.}
   \item{data.set}{Character string naming the marker data set to restrict.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper. The meaningful result is the side effect of restricting the marker data in the selected data set. Validation failures may return \code{NULL}.
}

\details{
The function modifies \code{data.set} in place, retaining only
individuals and markers that satisfy all applicable criteria.
A prior copy with \code{\link{st.copy.marker.data}} preserves the original
data set.

Markers are retained when their number of non-missing alleles does
not exceed \code{NoAll.MAX}, their missing-data proportion does not
exceed \code{MaMis.MAX}, and their expected heterozygosity is at
least \code{ExHet.MIN}. Individuals are retained when their
missing-data proportion does not exceed \code{InMis.MAX}.
Missingness thresholds are proportions, for example \code{0.1}
for 10 percent.

All criteria are evaluated using the statistics available at the
start of the call. Statistics are not recalculated between criteria
within the same call. Separate calls provide staged filtering.
For example, one call removes individuals with high missingness and
a later call filters markers using the updated statistics.

The linkage map and associated performance data are subset
along with the marker data. Removing markers invalidates
marker-dependent model results, which must be recalculated
before further use.
}
\note{
Restrictions are based on the marker statistics present at the start of each call. Repeating the function after statistics have been recalculated can therefore give a different result from applying several criteria simultaneously.
}

\examples{
data("v-tropmaize-vcf")
data("v-tropmaize-phe")

st.load.vcf.data(v.tropmaize.vcf)
st.load.performance.data(v.tropmaize.phe)

st.restrict.marker.data( NoAll.MAX = 2    )
st.restrict.marker.data( MaMis.MAX = 0.05 )
st.restrict.marker.data( ExHet.MIN = 0.1  )
st.restrict.marker.data( InMis.MAX = 0.1  )
}

%
