\name{evaluate.allele.freq2}
\alias{evaluate.allele.freq2}
\title{evaluate.allele.freq2()}

\description{
Calculates allele counts and relative frequencies at one or more loci
in a simulation population.
}

\usage{
evaluate.allele.freq2(PopName, Loci, missing=0)
}

\arguments{
   \item{PopName}{Character string naming the population to evaluate.}
   \item{Loci}{Character vector of locus names.}
   \item{missing}{Integer controlling missing-allele handling.
   With \code{0}, missing alleles are included.
   With \code{1}, missing allele copies are excluded from the counts
   and frequency calculation.}
}

\value{
A \code{data.frame} containing locus columns, \code{count}, and
\code{freq}.
}

\details{
The population and linkage map must already be defined.
Each element of \code{Loci} is interpreted as one locus name.
Space-separated names within a single string are not split.

For multiple loci, each row represents an allele at one locus;
the other locus columns contain zero placeholders.
Frequencies are calculated by dividing each count by the sum of
all returned counts. They therefore sum to one across the complete
table, not separately within each locus. To obtain frequencies
within a locus, evaluate that locus separately.
}

\examples{
map <- data.frame(
   chrom = 1:2,
   pos = c(0.01, 0.02),
   name = c("form", "color"),
   class = paste0("trait", 1:2)
)
define.genome(map)

init.population("P1", homozygote(1))
init.population("P2", homozygote(2))
cross("F1", "P1", "P2", 10)
cross("F2", "F1", "F1", 1000, self = 1)

# Allele frequencies at a single locus
evaluate.allele.freq2("F2", "form")

# Allele counts at both loci, with frequencies normalised
# across the complete table
evaluate.allele.freq2("F2", map$name)

# Evaluate each locus separately for within-locus frequencies
evaluate.allele.freq2("F2", "form")
evaluate.allele.freq2("F2", "color")
}

\seealso{
\code{\link{evaluate.allele.freq}}, an alternative implementation
that accepts space-separated locus names and differs in its
handling of missing alleles.
}