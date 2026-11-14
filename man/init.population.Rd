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
   \item{data}{Data frame with five numeric columns named \code{ind}, \code{chrom},
   \code{hom}, \code{pos}, and \code{all}, representing individual,
   chromosome, homologue, segment position, and allele.}
   \item{delete}{Logical. \code{TRUE} replaces an existing population with the same name
   after the input data have passed validation.}
}
\value{
\code{NULL}, invisibly. The function initializes the named simulation
population from the supplied chromosome-segment data.
}

\details{
The input data frame has exactly the five named numeric columns described above. Individual, chromosome and homologue numbers and allele codes are finite whole numbers; positions are finite and lie within the defined genome. Each row describes the start of a chromosome segment carrying the specified allele. The data must be ordered consistently with the simulation representation.

The number of input rows must be positive. Chromosome and homologue identifiers must be compatible with the currently defined simulation genome. The routine constructs the raw chromosome representation; marker genotypes and evaluated genetic values are derived data and are not implied merely by initialization.

Malformed input is rejected before an existing population is changed. Once
validation succeeds, \code{delete=TRUE} removes an existing population of
the same name before initialization. If creation fails while chromosome segments are being inserted, the incomplete new population is removed rather than retained.

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