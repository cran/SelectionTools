\name{sample.population}
\alias{sample.population}
\title{sample.population()}

\description{
Creates a new population by sampling individuals from an existing simulation population.
}

\usage{
sample.population(NameP1, NameP2, NoI=-1, rep=1)
}

\arguments{
   \item{NameP1}{Name of the sample population}
   \item{NameP2}{Name of the initial population}
   \item{NoI}{non-negative integer giving the number of individuals to choose}
   \item{rep}{Should sampling be with replacement?}
}
\value{
Returns invisible list of input arguments. The function is a compatibility wrapper and the meaningful result is the delegated population operation.
}

\details{
Each derived cache is complete for all individuals or absent throughout a
population. A mixed source cache causes all marker-genotype, genetic-value and
phenotype caches to be cleared silently in the source and result.

Sampling retains the source individuals; mixed source caches are cleared as described above. The requested sample is copied into the destination. Sampling can be with or without replacement according to \code{rep}. Raw chromosome data are always the authoritative copied representation; cached genotype and genetic-value categories are transferred only when complete in the source and successfully constructed for the entire destination.
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
sample.population("Sample", "Maize", NoI = 4, rep = 0)
population.size.get("Sample")
}

\seealso{
This function is a wrapper for \code{\link{population.sample}}.
}