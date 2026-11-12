\name{evaluate.allele2}
\alias{evaluate.allele2}
\title{evaluate.allele2()}

\description{
Calculates allele counts at one or more loci in a simulation population.
}

\usage{
evaluate.allele2(PopName, Loci, missing=0)
}

\arguments{
   \item{PopName}{Character string naming the population to evaluate.}
   \item{Loci}{Character vector of locus names. Each element specifies
   one locus.}
   \item{missing}{Integer controlling missing-allele handling.
   With \code{0}, missing alleles are included.
   With \code{1}, missing allele copies are excluded, while known
   copies are retained even if another allele of the same
   individual is missing.}
}

\value{
A \code{data.frame} containing locus columns and \code{count}.
Returns \code{NULL} invisibly when the request cannot be evaluated.
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
evaluate.allele2("Maize", loci[1:2])
}

\seealso{
\code{\link{evaluate.allele}} accepts space-separated locus names.
\code{\link{evaluate.allele.freq2}} also calculates relative frequencies.
}