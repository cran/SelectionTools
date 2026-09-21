\name{st.destroy.phase}
\alias{st.destroy.phase}
\title{Randomize the Stored Phase of a Marker Data Set}

\description{
Destroys known phase in a SelectionTools marker data set without changing any
unordered diploid genotype.  The function is mainly intended for validation of
statistical phasing on simulated data with known true phase.
}

\usage{
st.destroy.phase(data.set = "default")
}

\arguments{
  \item{data.set}{
    Character string naming an existing SelectionTools marker data set.
  }
}

\details{
For every observed heterozygous genotype, the two stored alleles are exchanged
with probability 0.5.  Homozygous and missing genotypes are unchanged.  Thus
the allele pair at every marker remains biologically identical, while the
correspondence of the two stored homologues across heterozygous markers is
randomized.

The function uses R's random-number generator.  Call \code{set.seed()} before
\code{st.destroy.phase()} when a reproducible randomized phase is required.
The marker data set is modified in place.
}

\value{
Invisibly returns a list with components \code{heterozygotes}, the number of
observed heterozygous genotypes examined, and \code{swapped}, the number whose
two stored alleles were exchanged.
}

\seealso{
\code{st.switch.error}, \code{st.phase}, \code{write.vcf}
}
