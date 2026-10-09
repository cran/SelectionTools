\name{st.get.simpop2}
\alias{st.get.simpop2}
\title{Return Simulation-Population Marker Data without a Linkage Map}
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
The transfer is restricted to diploid simulation data. Generated individual
identifiers use the same convention as \code{st.get.simpop()}: a population
P1 produces P1I1, P1I2, and so on. Population names and marker names
transferred from the simulation data must contain only ASCII letters and
digits. Invalid names are rejected rather than modified.

The two homologues returned by the simulation backend are chromosome-homologue
ordered. Consequently the resulting SelectionTools marker data set is marked
as phased.
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of copying simulation-population marker data without a linkage map into the named SelectionTools data set.
}
