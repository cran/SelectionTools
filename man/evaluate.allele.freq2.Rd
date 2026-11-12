\name{evaluate.allele.freq2}
\alias{evaluate.allele.freq2}
\title{evaluate.allele.freq2()}

\description{
Calculates allele counts and relative frequencies at one or more loci
in a simulation population.
}

\usage{
evaluate.allele.freq2(PopName, Loci, missing=0)
}

\arguments{
   \item{PopName}{Character string naming the population to evaluate.}
   \item{Loci}{Character vector of locus names.}
   \item{missing}{Integer controlling missing-allele handling.
   With \code{0}, missing alleles are included.
   With \code{1}, missing allele copies are excluded from the counts
   and frequency calculation.}
}

\value{
A \code{data.frame} containing locus columns, \code{count}, and
\code{freq}.
}

\details{
The population and linkage map must already be defined.
Each element of \code{Loci} is interpreted as one locus name.
Space-separated names within a single string are not split.

For multiple loci, each row represents an allele at one locus;
the other locus columns contain zero placeholders.
Frequencies are calculated by dividing each count by the sum of
all returned counts. They therefore sum to one across the complete
table, not separately within each locus. To obtain frequencies
within a locus, evaluate that locus separately.
}

\examples{
data("v-tropmaize-vcf", package = "SelectionTools")
marker.frame <- st.STvcf.to.dataframe(v.tropmaize.vcf)[1:12, c(1:6, 7:14)]
vcf <- st.dataframe.to.STvcf(marker.frame)
st.load.vcf.data(vcf, data.set = "maize.example")
reset.all()
st.set.sim.gp(data.set = "maize.example")
st.set.sim.mp(data.set = "maize.example")
st.set.sim.pp("Maize", data.set = "maize.example")
loci <- vcf$ID
evaluate.allele.freq2("Maize", loci[1])
evaluate.allele.freq2("Maize", loci[1:2])
}

\seealso{
\code{\link{evaluate.allele.freq}}, an alternative implementation
that accepts space-separated locus names and differs in its
handling of missing alleles.
}