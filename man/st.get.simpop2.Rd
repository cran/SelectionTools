\name{st.get.simpop2}
\alias{st.get.simpop2}
\title{Returns Simulation-Population Marker Data without a Linkage Map}
\description{
Imports diploid simulation-population marker genotypes into a SelectionTools
marker data set using the alternative simulation transfer path.
}
\usage{
st.get.simpop2(pop.name, data.set ="default")
}
\arguments{
\item{pop.name}{Character string naming the simulation population. For this
transfer the population name must contain 2 to 255 ASCII letters or digits.}
\item{data.set}{Character string naming the SelectionTools marker data set.}
}
\details{
The transfer supports only diploid simulation data. Stable individual names and full simulated pedigrees are retained, using the same rules as \code{st.get.simpop()}. Duplicated individual IDs are rejected before replacing the destination. Population names must contain only ASCII letters and digits. Marker names follow the same identifier rules as the
VCF interface: names must be non-empty, occupy fewer than 256 bytes, and
contain no spaces or ASCII control characters. Punctuation such as hyphens,
dots, and underscores is preserved. Invalid names are rejected rather than
modified.

The two homologues returned by the simulation backend are chromosome-homologue
ordered. Consequently the resulting SelectionTools marker data set is marked
as phased.
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of copying simulation-population marker data without a linkage map into the named SelectionTools data set.
}
