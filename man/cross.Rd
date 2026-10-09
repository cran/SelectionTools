\name{cross}
\alias{cross}
\title{cross()}

\description{
Generates a diploid progeny population by crossing individuals drawn from two parental populations.
}

\usage{
cross(NamePg, NameP1, NameP2, NoPg = -1, self = 0, mode = 1,
      crossing.scheme = 0)
}

\arguments{
\item{NamePg}{Character string naming the progeny population.}
   \item{NameP1}{Character string naming parental population 1.}
   \item{NameP2}{Character string naming parental population 2.}
   \item{NoPg}{Number of progeny to generate. A positive value creates that many individuals; \code{-1} uses the allocation of an existing progeny population.}
  \item{self}{Probability of explicitly selfing a selected parent,
  between 0 and 1. With probability \code{self/2}, the selected
  individual from parental population 1 is selfed; with the same
  probability, the selected individual from parental population 2
  is selfed. Otherwise, the two selected individuals are crossed.
  Thus, \code{self=1} produces only selfed progeny.}
  \item{mode}{Integer controlling whether the same individual may
  be selected as both parents before explicit selfing is applied.
  With \code{mode=0}, the second parent is resampled if both selected
  parents are the same individual and the parental population contains
  more than one individual. With \code{mode=1}, this restriction is
  not applied.}
   \item{crossing.scheme}{Integer 0 to 3 selecting random or deterministic parent choice in the two parental populations.}
}
\value{
Invisibly returns a list containing the arguments passed to the
compiled crossing routine. The function is primarily called for
its side effect of creating the progeny population.
}

\details{ The simulation genome must be diploid. The progeny population
must have a name different from both parental populations. 
Both parental population names may refer to the same population.
Using the same population for both parents does not, by itself,
ensure selfing. Existing contents under the progeny name are replaced when a new progeny
population is constructed.

For positive \code{NoPg}, exactly that many progeny individuals are
generated. The special value \code{-1} uses the already allocated size
of an existing progeny population; it is therefore meaningful only when
such a destination population already exists with positive allocation.

The four crossing schemes control whether the parent chosen from each
parental population is random or deterministic. Deterministic selection
cycles through the active individuals of the corresponding parent
population. For progeny number \eqn{n} (counted from 1), the one-based
position of the deterministically selected parent in parental population
\eqn{P} of size \eqn{S_P} is
\deqn{i_P(n) = ((n-1) \bmod S_P) + 1.}
Thus deterministic selection cycles through all active individuals and
starts again with the first individual after the end of the parental
population is reached. The rule is applied separately to each parental
population, so the two populations may have different sizes.
\code{self} controls selfing, and \code{mode} controls whether the mating
logic applies self-incompatibility.

Each progeny receives newly simulated chromosome segments produced by
meiosis. Cached marker genotypes and evaluated genetic values are not
inherited as authoritative derived data; they can be generated from the
new chromosomes when needed. If chromosome construction fails, the
incomplete progeny population is removed rather than left partially
generated.  }

\note{ Crossing scheme 0 selects both parents randomly; scheme 1 selects
the first parent randomly and the second deterministically; scheme 2
does the reverse; scheme 3 selects both deterministically.  }

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
}