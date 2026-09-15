\name{st.set.sim.ef}
\alias{st.set.sim.ef}
\title{Transfer Marker Effects to the Simulation Backend}

\description{
Transfers estimated marker effects from a SelectionTools marker data set
directly to the active simulation backend.
}

\usage{
st.set.sim.ef(name="effect", data.set="default")
}

\arguments{
\item{name}{Character string naming the effect in the simulation backend.}
\item{data.set}{Character string naming the SelectionTools marker data set
that provides the estimated marker effects.}
}

\details{
The function transfers marker effects directly between the genomic-selection
and simulation data structures; no intermediate effect file is written.

With the default allele codes used by the genomic-selection routines,
\code{gs.set.allele.codes()} assigns the additive dosage codes 0, 1, and 2
to the genotypes with zero, one, and two copies of the effect allele,
respectively. For marker \eqn{j} with estimated effect \eqn{u_j}, the
marker contribution in the genomic-selection representation is therefore
\eqn{0}, \eqn{u_j}, or \eqn{2u_j}. The default codes for genotypes
containing a missing allele are 0.5, 1, and 1.5.

The simulation backend stores an additive effect per allele. For each
transferred marker, the effect allele receives \eqn{+u_j/2} and the
complementary allele receives \eqn{-u_j/2}. Before adjustment of the
general mean, the three diploid genotype contributions are consequently
\eqn{-u_j}, \eqn{0}, and \eqn{+u_j}. This representation is equivalent to
changing the marker dosage from \eqn{z_{ij}} to
\deqn{z^*_{ij}=z_{ij}-1.}{z*_[ij] = z_[ij] - 1.}

For the genomic-selection representation, the estimated genetic value is
\deqn{
\hat{g}_i=\hat{\mu}_{012}+\sum_{j=1}^{m}z_{ij}\hat{u}_j.
}{
g_[i] = mu_[012] + sum(z_[ij] u_[j]).
}
Using \eqn{z_{ij}=z^*_{ij}+1} gives
\deqn{
\hat{g}_i=
\left(\hat{\mu}_{012}+\sum_{j=1}^{m}\hat{u}_j\right)
+\sum_{j=1}^{m}z^*_{ij}\hat{u}_j.
}{
g_[i] = (mu_[012] + sum(u_[j])) + sum(z*_[ij] u_[j]).
}
Thus, the general mean stored for the simulation effect is
\deqn{
\hat{\mu}_{\mathrm{sim}}=
\hat{\mu}_{012}+\sum_{j=1}^{m}\hat{u}_j.
}{
mu_[sim] = mu_[012] + sum(u_[j]).
}
This is the coding-shift correction described by Stranden and
Christensen (2011), Equation (2). In their general notation the change in the
estimated general mean is the coding-shift vector transposed times the
marker-effect vector. For the shift from 0/1/2 to -1/0/1, the coding-shift
vector consists of ones and the correction is therefore the sum of marker
effects.

This shift is distinct from allele-frequency centering. If a marker column is
centered as \eqn{z^c_{ij}=z_{ij}-2p_j}, the corresponding change of the
general mean is \eqn{2\sum_j p_j u_j}. The simulation transfer uses
\eqn{z^*_{ij}=z_{ij}-1}, not allele-frequency centering. The two shifts are
the same for a marker only when \eqn{p_j=0.5}.

Only markers that are present in the active simulation map are transferred.
Consequently, the sum in the intercept correction is taken over transferred
markers only. Markers that are absent from the simulation map are ignored
with a warning.

The exact equivalence described above applies to the default allele codes.
\code{gs.set.allele.codes()} permits power users to define other codes; an
arbitrary custom coding is not in general represented by the fixed
\eqn{-u_j/2}, \eqn{+u_j/2} simulation allele effects.

Marker effects are optional. If \code{data.set} contains no estimated marker
effects, the function returns successfully without changing the simulation
effects. A partially inconsistent marker-effect state is treated as an error.

When marker effects are present, a diploid simulation genome and a linkage
map must already be installed. Marker effects are connected to simulation
loci by marker name. The data set supplying the effects therefore need not be
the data set that originally supplied the simulation map.

If an effect with the requested \code{name} is already loaded, the operation
fails. After a new effect is installed successfully, stored population
genetic and phenotypic values are removed because their effect dimension is
no longer current.
}

\value{
An invisible \code{list} returned by the underlying \code{.C} routine. Component \code{r} is the status value used by the simulation bridge; the function is primarily called for its side effect of transferring the requested data to the simulation backend.
}

\references{
Stranden I, Christensen OF (2011) Allele coding in genomic evaluation.
Genetics Selection Evolution 43:25. See Equation (2).
}

\seealso{
\code{gs.set.allele.codes}, \code{st.set.sim.gp}, \code{st.set.sim.mp},
\code{st.set.sim.pp}, \code{st.set.simpop}
}
