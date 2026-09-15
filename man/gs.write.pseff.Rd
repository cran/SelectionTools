\name{gs.write.pseff}
\alias{gs.write.pseff}
\title{Write Marker Effects for the Simulation Backend}

\description{
Writes estimated marker effects in the effect-file format used by the
simulation backend.
}

\usage{
gs.write.pseff(e, beta0 = NA, file)
}

\arguments{
\item{e}{A marker-effect table, typically returned by
\code{gs.return.effects()}. The first row contains the general mean and the
remaining rows contain marker names, effect alleles, complementary alleles,
and estimated marker effects.}
\item{beta0}{Optional general mean on the genomic-selection scale. If
\code{NA}, the general mean in the first row of \code{e} is used.}
\item{file}{Character string naming the effect file to be written. This argument is required because the function writes a persistent file.}
}

\details{
The genomic-selection and simulation parts of SelectionTools use different
origins for additive marker coding.

With the default allele codes used by the genomic-selection routines,
\code{gs.set.allele.codes()} assigns the additive dosage codes 0, 1, and 2
to the genotypes with zero, one, and two copies of the effect allele,
respectively. For marker \eqn{j} with estimated effect \eqn{u_j}, the
marker contribution is therefore \eqn{0}, \eqn{u_j}, or \eqn{2u_j}.

The simulation effect-file representation stores an additive effect per
allele. The effect allele is written with effect \eqn{+u_j/2} and the
complementary allele with effect \eqn{-u_j/2}. The corresponding diploid
marker contributions are therefore \eqn{-u_j}, \eqn{0}, and \eqn{+u_j}.
This is equivalent to changing the marker dosage from \eqn{z_{ij}} to
\deqn{z^*_{ij}=z_{ij}-1.}{z*_[ij] = z_[ij] - 1.}

On the genomic-selection scale,
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
Consequently, the first line written to the effect file is
\deqn{
\hat{\mu}_{\mathrm{sim}}=
\hat{\mu}_{012}+\sum_{j=1}^{m}\hat{u}_j.
}{
mu_[sim] = mu_[012] + sum(u_[j]).
}
If \code{beta0} is supplied, it replaces \eqn{\hat{\mu}_{012}} before this
coding correction is applied.

This is the coding-shift correction described by Stranden and
Christensen (2011), Equation (2). The shift is distinct from
allele-frequency centering. For coding
\eqn{z^c_{ij}=z_{ij}-2p_j}, the corresponding mean shift is
\eqn{2\sum_j p_j u_j}; the simulation effect-file coding instead uses
\eqn{z^*_{ij}=z_{ij}-1}.

The exact equivalence described above applies to the default allele codes.
\code{gs.set.allele.codes()} permits power users to define other codes; an
arbitrary custom coding is not in general represented by the fixed
\eqn{-u_j/2}, \eqn{+u_j/2} simulation allele effects.
}

\value{
\code{NULL}, returned invisibly by the file-writing operation. The function is called for its side effect of writing marker effects for the simulation backend.
}

\references{
Stranden I, Christensen OF (2011) Allele coding in genomic
evaluation. Genetics Selection Evolution 43:25. See Equation (2).
}

\seealso{
\code{gs.return.effects}, \code{gs.set.allele.codes},
\code{st.set.sim.ef}, \code{load.effmap}
}
