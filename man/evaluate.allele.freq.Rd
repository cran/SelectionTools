\name{evaluate.allele.freq}
\alias{evaluate.allele.freq}
\title{Calculate or return allele frequencies for selected loci}

\description{
Calculate or return allele frequencies for selected loci.
}

\usage{
evaluate.allele.freq(PopName, Loci, missing=0)
}

\arguments{
   \item{PopName}{Name of the population.}
    \item{Loci}{A character string containing one locus name or
    multiple locus names separated by spaces. Character vectors
    of length greater than one are not supported.}
    
    \item{missing}{Integer controlling missing-allele handling.
    With \code{0}, missing alleles are included.
    With \code{1}, missing alleles are excluded.}
}
\value{
A \code{data.frame} containing the allele counts and relative frequencies.
}

\examples{
map <- data.frame(chrom = 1,
                  pos = 0.01,
                  name = "form")
define.genome(map)

init.population("P1", homozygote(1))
init.population("P2", homozygote(2))

# Example: Generate 100 F1 hybrids
cross("F1","P1", "P2", 100)
evaluate.genotype("P1","form")

# Example: Repeated selfing
cross("F2", "F1", "F1", 1000, self = 1)
cross("F3", "F2", "F2", 1000, self = 1)
cross("F4", "F3", "F3", 1000, self = 1)

evaluate.genotype("F4", "form")
evaluate.allele.freq("F4","form")
}

\seealso{
For vector based input of loci, see \code{\link{evaluate.allele.freq2}}.
}