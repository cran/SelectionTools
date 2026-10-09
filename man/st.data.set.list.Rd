% copyright (C) 2001 Frisch/Maurer
%
\name{st.data.set.list}
\alias{st.data.set.list}
\alias{print.st.data.set.list}
\title{List SelectionTools Marker Data Sets}
\description{
Returns a compact summary of the marker data sets currently registered in
SelectionTools.
}
\usage{
st.data.set.list(digits=getOption("digits"))

\method{print}{st.data.set.list}(x, digits=attr(x,"digits"), complete=FALSE, ...)
}
\arguments{
\item{digits}{Integer between 1 and 22 giving the number of significant digits
used when the returned summary is printed. The stored numeric values are not
rounded.}
\item{x}{An object returned by \code{st.data.set.list()}.}
\item{complete}{Logical. If \code{FALSE} (the default), print only the compact overview columns \code{data.set}, \code{no.markers}, \code{no.individuals}, \code{map}, \code{ram.core}, and \code{ram.work}. If \code{TRUE}, print all columns.}
\item{...}{Additional arguments passed to \code{print.data.frame()}.}
}
\details{
Only non-empty data sets are returned. A data set is omitted when both its
number of markers and number of individuals are zero. Thus the built-in
\code{"default"} data set is not shown while it is empty. Non-empty
user-defined data sets are returned in creation order.

\code{heterozygosity} is the fraction of complete diploid genotype calls that
are heterozygous. \code{missing.markers} is the fraction of
marker-by-individual genotype calls with at least one missing allele. Both are
between zero and one and are \code{NA} when they cannot be calculated.

\code{phased} is \code{TRUE} when the stored \code{All1}/\code{All2} order is
known to be phased. VCF input is marked phased only when genotype separators
are consistently \code{|}; STvcf import/export retains the same data-set phase
state. Simulation populations returned by \code{st.get.simpop()} or
\code{st.get.simpop2()} are phased because the simulation backend returns the
two chromosome homologues separately. A complete successful \code{st.phase()}
commit is marked phased when no heterozygous phase remains unresolved, while
\code{st.destroy.phase()} clears the state. Operations that only subset
individuals or markers, replace a linkage map, add phenotypic values, copy
marker data, or replace markers by derived haplotype blocks retain the state.

\code{map} indicates that a linkage map is present and
\code{no.chromosomes} gives the number of chromosomes represented by that map.
\code{effects}, \code{phenotypic.value}, and \code{genotypic.value} indicate
the presence of estimated effects, phenotypic values, and genotypic values,
respectively.

\code{ram.core} and \code{ram.work} are estimates in MiB (1024^2 bytes).
\code{ram.core} is deliberately restricted to the marker data required to
retain the data set: individual and marker names, the two allele arrays, and
the linkage-map storage. \code{ram.work} contains all other dynamically owned
state, including marker statistics, phenotypic and genotypic values, effects
and model results, haplotype/cross objects, chromosome-specific recombination
and LD matrices, and reusable calculation buffers. Allocator overhead is not
included. \code{st.compact.data.set()} releases the work portion while
retaining the core.

The returned object always contains all columns. Its print method shows a compact
overview by default: \code{data.set}, \code{no.markers},
\code{no.individuals}, \code{map}, \code{ram.core}, and \code{ram.work}.
Use \code{print(x,complete=TRUE)} to display the complete summary.

The \code{digits} argument affects printing only. Thus, for example,
\code{x <- st.data.set.list(digits=3)} stores full-precision values but prints
the numeric columns with three significant digits. Printing can also be
overridden later with \code{print(x,digits=5)}.
}
\value{
An object inheriting from \code{data.frame}, with one row per non-empty
registered data set and the following columns. If all registered data sets are
empty, a zero-row data frame with these columns is returned:
\itemize{
  \item \code{data.set}: data-set name.
  \item \code{no.markers}: number of markers.
  \item \code{no.individuals}: number of individuals.
  \item \code{heterozygosity}: heterozygosity among complete genotype calls.
  \item \code{missing.markers}: fraction of missing marker-genotype calls.
  \item \code{phased}: whether allele order is known to be phased.
  \item \code{map}: whether a linkage map is present.
  \item \code{no.chromosomes}: number of chromosomes represented in the map.
  \item \code{effects}: whether estimated effects are present.
  \item \code{phenotypic.value}: whether phenotypic values are present.
  \item \code{genotypic.value}: whether genotypic values are present.
  \item \code{ram.core}: estimated marker-data and linkage-map RAM in MiB.
  \item \code{ram.work}: estimated calculated/result/workspace RAM in MiB.
}
The \code{digits} setting is stored as an attribute and changes printing only.
}
