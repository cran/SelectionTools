\name{st.set.sim.pp}
\alias{st.set.sim.pp}
\title{Transfer a Marker Population to the Simulation Backend}

\description{
Creates a simulation population from the diploid marker matrix of a
SelectionTools data set.
}

\usage{
st.set.sim.pp(pop.name="default", data.set="default")
}

\arguments{
\item{pop.name}{Character string naming the simulation population to create.
The name must contain 2 to 255 ASCII letters or digits.}
\item{data.set}{Character string naming the SelectionTools marker data set
supplying individuals and alleles.}
}

\details{
A diploid simulation genome and linkage map must already be installed. The
source marker data set does not need to contain the linkage map that was used
to initialize the simulation. Instead, the active simulation loci are matched
to rows of the source marker matrix by marker name.

Every marker in the active simulation linkage map must occur in the source
marker matrix. The source data set may contain additional markers; these are
ignored during population transfer. Consequently, populations from several
compatible SelectionTools data sets can be added to one initialized simulation
without reloading genome parameters or the linkage map.

Imported individual identifiers and full stored pedigrees are copied into owned simulation strings. If no pedigree is attached, the imported name is the founder pedigree. Removing or changing the source marker dataset does not change those strings. No fitted symbolic-A prediction state is transferred.

For each imported individual, allele values are placed at the corresponding
simulation loci. SelectionTools biological allele codes 1 through 254 are
transferred unchanged, while the internal SelectionTools missing representation
is converted to the simulation missing-allele indicator. After construction,
the imported population is checked marker by marker against the source marker matrix. If
verification fails, the newly created population is removed.

Compatibility of biological coding across independently supplied data sets is
the responsibility of the user. In particular, matching marker names are
assumed to identify the same loci and allele coding. The requested
\code{pop.name} must not already exist in the simulation backend.
}

\value{
An invisible \code{list} returned by the compiled routine. Component \code{r} is the status value used by the simulation bridge; the function is primarily called for its side effect of transferring the requested data to the simulation backend.
}

\seealso{
\code{st.set.sim.gp}, \code{st.set.sim.mp}, \code{st.set.sim.ef},
\code{st.set.simpop}
}
