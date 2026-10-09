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
   \item{delete}{Logical. If true, remove an existing population with the same name before initialization.}
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
# Define a population explicitly
set.genome.par(no.chrom = 1, no.hom = 2, chrom.len = 1)
define.map("form trait1 1/0.5")

# Individual 1: genotype 1/1
# Individual 2: genotype 1/2
pop <- data.frame(
   ind   = c(1, 1, 2, 2),
   chrom = c(1, 1, 1, 1),
   hom   = c(1, 2, 1, 2),
   pos   = c(0, 0, 0, 0),
   all   = c(1, 1, 1, 2)
)

init.population("P1", pop)
evaluate.genotype("P1", "form")

# Define a genome from a data frame and create homozygous individuals
reset.all()

map <- data.frame(
   chrom = 1,
   pos = 0.5,
   name = "form"
)
define.genome(map)

# Two individuals, both homozygous for allele 1
init.population("P2", homozygote(1, NoInd = 2))
evaluate.genotype("P2", "form")
}