\name{append.population}
\alias{append.population}
\title{append.population()}

\description{
Appends a copy of one population to another.
}

\usage{
append.population(NameP1, NameP2)
}

\arguments{
   \item{NameP1}{Character string naming the destination population.
   If it does not exist, it is created.}
   \item{NameP2}{Character string naming the source population.
   This population is not modified.}
}

\value{
Returns NameP1 and NameP2. Usually called for side effects.
}

\details{
The individuals of \code{NameP2} are copied to the end of
\code{NameP1}. Existing individuals in \code{NameP1} are retained,
and \code{NameP2} is not modified.

The source population must exist and contain at least one individual.
The source and destination population names must differ.
}

\examples{
# Define the genome and linkage map
map <- data.frame(
   chrom = 1,
   pos = 0.5,
   name = "form",
   class = "trait1"
)
define.genome(map)

# Create two populations with different alleles
init.population("P1", homozygote(1, NoInd = 3))
init.population("P2", homozygote(2, NoInd = 2))

# Combine copies of both populations in a new population
append.population("Combined", "P1")
append.population("Combined", "P2")

# Combined contains three individuals with genotype 1/1
# and two individuals with genotype 2/2.
# The original populations P1 and P2 are retained.
evaluate.genotype("Combined", "form")
population.list()
}

\seealso{
\code{\link{population.append}}, which is called by this wrapper.
}