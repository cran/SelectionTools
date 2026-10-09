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
map <- data.frame(
   chrom = 1:2,
   pos = c(0.01, 0.02),
   name = c("form", "color")
)
define.genome(map)
init.population("P1", homozygote(1, NoInd = 2))

evaluate.allele("P1", "form color")
}

\seealso{
\code{\link{evaluate.allele2}} accepts a character vector of locus
names and excludes missing allele copies without first removing
loci lacking completely observed genotypes.
\code{\link{evaluate.allele.freq}} also calculates relative frequencies.
}