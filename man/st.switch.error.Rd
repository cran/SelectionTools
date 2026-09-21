\name{st.switch.error}
\alias{st.switch.error}
\title{Calculate Switch Error against Known Phase}

\description{
Compares a phased SelectionTools marker data set with a data set containing
known true phase and calculates the local switch-error rate chromosome by
chromosome.
}

\usage{
st.switch.error(true.data.set, data.set = "default")
}

\arguments{
  \item{true.data.set}{
    Character string naming the SelectionTools marker data set containing the
    true ordered homologues.
  }
  \item{data.set}{
    Character string naming the phased data set to evaluate.
  }
}

\details{
Individuals and markers are matched by their unique names, so the two data
sets need not use the same storage order.  They must contain the same numbers
and names of individuals and markers and complete linkage maps.

At each observed heterozygote in the true data, the candidate phase is recorded
as having the same or opposite orientation relative to the truth.  Two
consecutive scorable heterozygotes on the same chromosome define one switch
link.  A switch error is counted when their relative orientations differ.
The orientation state is reset at chromosome boundaries, so an arbitrary
exchange of the two haplotype labels for an entire chromosome has no effect on
the result.  Missing candidate genotypes and incompatible candidate genotypes
are not used to form links and are reported separately.

The current implementation is for biallelic markers, matching the current
\code{st.phase()} implementation.  Allele codes are reduced to marker-specific
ranks before comparison, which permits direct comparison with data exported by
\code{write.vcf()}, phased by Beagle, and re-imported by \code{read.vcf()}.

The switch-error rate is the number of switch errors divided by the number of
scored links.  This is the standard local switch-error definition used in
statistical-phasing evaluations; see Browning and Browning (2007) and Browning
and Browning (2022).
}

\value{
A data frame with one column, \code{result}.  The row names are
\code{errors}, \code{links}, \code{switch.error},
\code{switch.error.percent}, \code{heterozygotes},
\code{scored.heterozygotes}, \code{coverage},
\code{missing.heterozygotes}, and \code{genotype.mismatches}.
\code{switch.error} and \code{coverage} are proportions; the former is
\code{NA} when no switch link is available and the latter is \code{NA} when
the true data contain no observed heterozygotes.  The one-column layout is
intended to make the result directly readable when printed.
}

\references{
Browning SR, Browning BL (2007) Rapid and accurate haplotype phasing and
missing-data inference for whole-genome association studies by use of localized
haplotype clustering. \emph{American Journal of Human Genetics} 81:1084--1097.
\doi{10.1086/521987}.

Browning BL, Browning SR (2022) Genotype error biases trio-based estimates of
haplotype phase accuracy. \emph{American Journal of Human Genetics}
109:1016--1025. \doi{10.1016/j.ajhg.2022.04.019}.
}

\seealso{
\code{st.destroy.phase}, \code{st.phase}, \code{read.vcf}, \code{write.vcf}
}
