\name{evaluate.allele}
\alias{evaluate.allele}
\title{evaluate.allele()}

\description{
Calculates allele counts at one or more loci in a simulation population.
}

\usage{
evaluate.allele(PopName, Loci, missing=0)
}

\arguments{
   \item{PopName}{Character string naming the population to evaluate.}
   \item{Loci}{A character string containing locus names separated
   by spaces.}
   \item{missing}{Integer controlling missing-allele handling.
   With \code{0}, missing alleles are included.
   With \code{1}, loci without any completely observed genotype
   are removed, and missing allele copies at the remaining loci
   are excluded from the counts.}
}

\value{
A \code{data.frame} containing locus columns and \code{count}.
Returns \code{NULL} invisibly when no result is available.
}

\details{
The population and linkage map must already be defined.
For multiple loci, each row represents an allele at one locus.
The other locus columns contain zero placeholders.
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
evaluate.allele("Maize", loci[1])
evaluate.allele("Maize", paste(loci[1:2], collapse = " "))
}

\seealso{
\code{\link{evaluate.allele2}} accepts a character vector of locus
names and excludes missing allele copies without first removing
loci lacking completely observed genotypes.
\code{\link{evaluate.allele.freq}} also calculates relative frequencies.
}