\name{cross}
\alias{cross}
\alias{population.cross}
\title{cross()}

\description{
Generates a diploid progeny population by crossing individuals drawn from two parental populations.
}

\usage{
cross(NamePg, NameP1, NameP2, NoPg = -1, self = 0, mode = 1,
      crossing.scheme = 0)

population.cross(name, pool1, pool2, n)
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


\item{name}{Name of the offspring simulation population.}
\item{pool1,pool2}{Names of two distinct, nonempty parental simulation populations. Neither may be named \code{"name"}.}
\item{n}{Positive integer number of crosses.}

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
data("v-tropmaize-vcf", package = "SelectionTools")
marker.frame <- st.STvcf.to.dataframe(v.tropmaize.vcf)[1:12, c(1:6, 7:14)]
vcf <- st.dataframe.to.STvcf(marker.frame)
st.load.vcf.data(vcf, data.set = "maize.example")
reset.all()
st.set.sim.gp(data.set = "maize.example")
st.set.sim.mp(data.set = "maize.example")
st.set.sim.pp("Maize", data.set = "maize.example")
loci <- vcf$ID
cross("Offspring", "Maize", "Maize", 4)
population.size.get("Offspring")
}
\section{Crosses with a returned crossing design}{
Creates random crosses between two simulation populations and returns the selected parental identifiers at creation time.

Uses the existing diploid random-cross routine and its simulation random-number generator. Each offspring owns its full simulated pedigree independently of the returned design. The design records the immediate crossing keys; it is not reconstructed by parsing pedigree text. No genetic or phenotypic values are transferred. An existing offspring population is replaced, as by \code{cross}.

A data frame with \code{name} and two columns named after the parental populations. Rows follow offspring order and contain its stable identifier and the two selected parental identifiers.

}
