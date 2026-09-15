\name{st.get.simpop}
\alias{st.get.simpop}
\title{Return Simulation-Population Data in SelectionTools Format}
\description{
Imports a diploid simulation population into a SelectionTools marker data set.
}
\usage{
st.get.simpop(pop.name, data.set ="default")
}
\arguments{
\item{pop.name}{Character string naming the simulation population to import.
For this transfer the population name must contain 2 to 255 ASCII letters or
digits.}
\item{data.set}{Character string naming the SelectionTools marker data set to
create or replace.}
}
\details{
The simulation population is converted to the two-allele marker matrix
representation used by SelectionTools. The simulation genome must be diploid.
Marker identities and loci are obtained from the active simulation map.
Marker names transferred from that map must contain only ASCII letters and
digits; invalid names are rejected rather than modified. Simulation missing
alleles are converted to the SelectionTools internal missing representation,
and marker statistics are calculated for the imported matrix. Non-missing
simulation allele codes must lie between 1 and 254 to be transferred into a
SelectionTools marker data set; wider simulation allele codes are rejected
rather than narrowed or wrapped.

The transfer creates new marker-data individual identifiers rather than
recovering names that may have existed before the population entered
simulation. For a population named P1, the generated identifiers are P1I1,
P1I2, and so on. The separator I is alphanumeric so the resulting names obey
the marker-data/STvcf convention that individual identifiers contain only
letters and digits. Population names containing dots, punctuation,
whitespace, or other non-alphanumeric characters are rejected by this
transfer rather than being silently sanitized.
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of copying a simulation population into the named SelectionTools marker data set.
}
