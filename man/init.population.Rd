\name{init.population}
\alias{init.population}
\title{init.population()}

\description{
Creates a simulation population from the internal chromosome-segment data frame representation.
}

\usage{
init.population(name, data, delete=TRUE)
}

\arguments{
\item{name}{Character string naming the simulation population to initialize.}
   \item{data}{Data frame with five columns representing individual, chromosome, homologue, segment position, and allele.}
   \item{delete}{Logical. TRUE removes an existing population with the same name before initialization.}
}
\value{
The function initialises the named simulation population from the supplied data in the backend. 
Returns list of name, total number of alleles and input data columns.
The function is primarily called for its side effect of creating
the simulation population.
}

\details{
The input data frame must have exactly five columns named or interpretable as individual number, chromosome number, homologue number, position, and allele. Each row describes the start of a chromosome segment carrying the specified allele. The data must be ordered consistently with the simulation representation.

The number of input rows must be positive. Chromosome and homologue identifiers must be compatible with the currently defined simulation genome. The routine constructs the raw chromosome representation; marker genotypes and evaluated genetic values are derived data and are not implied merely by initialization.

When \code{delete=TRUE}, an existing population of the same name is removed before initialization. If creation fails while chromosome segments are being inserted, the incomplete new population is removed rather than retained.

Uses the currently defined genome and linkage map. Both must be defined before calling this function.
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
observed.population <- get.population("Maize")
init.population("Imported", observed.population)
evaluate.genotype("Imported", loci[1])
}