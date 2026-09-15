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
   \item{self}{Numeric selfing rate passed to the crossing routine.}
   \item{mode}{Integer mating-mode flag controlling self-incompatibility behavior.}
   \item{crossing.scheme}{Integer 0 to 3 selecting random or deterministic parent choice in the two parental populations.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of creating the requested crosses.
}

\details{ The simulation genome must be diploid. The progeny population
must have a name different from both parental populations. Existing
contents under the progeny name are replaced when a new progeny
population is constructed.

For positive \code{NoPg}, exactly that many progeny individuals are
generated. The special value \code{-1} uses the already allocated size
of an existing progeny population; it is therefore meaningful only when
such a destination population already exists with positive allocation.

The four crossing schemes control whether the parent chosen from each
parental population is random or deterministic. Deterministic selection
cycles through the active individuals of the corresponding parent
population. \code{self} controls selfing, and \code{mode} controls
whether the mating logic applies self-incompatibility.

Each progeny receives newly simulated chromosome segments produced by
meiosis. Cached marker genotypes and evaluated genetic values are not
inherited as authoritative derived data; they can be generated from the
new chromosomes when needed. If chromosome construction fails, the
incomplete progeny population is removed rather than left partially
generated.  }

\note{ Crossing scheme 0 selects both parents randomly; scheme 1 selects
the first parent randomly and the second deterministically; scheme 2
does the reverse; scheme 3 selects both deterministically.  }
