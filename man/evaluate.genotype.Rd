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
# Simple functionalty:
map <- data.frame(chrom = 1,
                  pos = 0.0,
                  name = "form",
                  class = "trait1")
define.genome(map)

init.population("P1", homozygote(1))
init.population("P2", homozygote(2))

cross("F1","P1", "P2", 100)
cross("BC1","F1", "P2", 100)

evaluate.genotype("BC1","form")

# Complex functionality: Two linked loci
reset.all()

map <- data.frame(
  chrom = c(1, 1),
  pos = c(0.2, 0.8),
  name = c("A", "B")
)
define.genome(map, length.vec = 1)

# Individual 1 carries haplotypes 1-1 and 2-2.
# Individual 2 carries haplotypes 1-2 and 2-1.
pop <- data.frame(
  ind   = c(1, 1, 1, 1, 2, 2, 2, 2),
  chrom = c(1, 1, 1, 1, 1, 1, 1, 1),
  hom   = c(1, 1, 2, 2, 1, 1, 2, 2),
  pos = c(0.2, 0.8, 0.2, 0.8, 0.2, 0.8, 0.2, 0.8),  
  all   = c(1, 2, 1, 2, 1, 2, 2, 1)
)
init.population("PhaseExample", pop)

# Preserve phase: the two individuals form separate categories,
# each with count 1 and frequency 0.5.
evaluate.genotype("PhaseExample", c("A", "B"), mode = 0)

# Ignore phase: both individuals are heterozygous at both loci,
# giving one category with count 2 and frequency 1.
evaluate.genotype("PhaseExample", c("A", "B"), mode = 2)

# Loci can also be supplied as a space-separated string.
evaluate.genotype("PhaseExample", "A B", mode = 2)

# Retrieve one observed phased combination.
evaluate.genotype("PhaseExample", c("A", "B"),
                  alleles = c(1, 2, 1, 2), mode = 0)
}