\name{evaluate.genotype}
\alias{evaluate.genotype}
\title{evaluate.genotype()}

\description{
  Calculates genotype frequencies
}

\usage{
evaluate.genotype(PopName, Loci, alleles=0, mode=2, bg=0)
}

\arguments{
   \item{PopName}{Name of the population to be evaluated}
   \item{Loci}{A List of Loci}
   \item{alleles}{A list of different Alleles}
   \item{mode}{Considering gametic phase or not}
   \item{bg}{Integer controlling missing-allele handling.
   With \code{0}, genotype combinations containing missing alleles
   are included. With \code{1}, combinations containing the missing
   allele indicator are excluded and frequencies are recalculated
   among the remaining individuals.}
}

\value{
A \code{data.frame} containing the genotype combinations and the count/frequency information returned by the compiled genotype evaluator. 
If a specific allele combination is requested the result can be reduced to the matching row. 
Returns \code{NULL} invisibly if no result is available.
}

\details{
The population and linkage map must already be defined.
Counts represent numbers of individuals,
and frequencies are calculated relative to the total number of
included individuals. Names must only consist of letters and numbers (no spaces, dots, etc.).
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
evaluate.genotype("Maize", loci[1])
evaluate.genotype("Maize", paste(loci[1:2], collapse = " "), bg = 1)
}