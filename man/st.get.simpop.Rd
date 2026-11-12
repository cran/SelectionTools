\name{st.get.simpop}
\alias{st.get.simpop}
\title{Returns Simulation-Population Data in SelectionTools Format}
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
Marker names transferred from that map follow the same identifier rules as
the VCF interface: names must be non-empty, occupy fewer than 256 bytes, and
contain no spaces or ASCII control characters. Punctuation such as hyphens,
dots, and underscores is preserved; invalid names are rejected rather than
modified. Simulation missing
alleles are converted to the SelectionTools internal missing representation,
and marker statistics are calculated for the imported matrix. Non-missing
simulation allele codes must lie between 1 and 254 to be transferred into a
SelectionTools marker data set; wider simulation allele codes are rejected
rather than narrowed or wrapped.

The transfer retains each individual's stable identifier and complete stored simulated pedigree. Imported founders retain their imported identifiers. New cross, DH and SSD individuals retain their creation names, including their full pedigree expressions. Copies and population renames keep these identities. A simulation population containing repeated copies with the same identifier cannot be exported to a marker dataset; the destination is left intact. \code{\link{population.pedigree}} reports names and pedigrees in either representation.


The two homologues returned by the simulation backend are chromosome-homologue
ordered. Consequently the resulting SelectionTools marker data set is marked
as phased.
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of copying a simulation population into the named SelectionTools marker data set.
}
