\name{dh}
\alias{dh}
\title{dh()}
\description{
Generates a diploid population of doubled-haploid lines from a parental population.
}
\usage{
dh(NamePg, NameP1, NoPg=-1)
}

\arguments{
\item{NamePg}{Character string naming the doubled-haploid progeny population.}
   \item{NameP1}{Character string naming the parental population.}
   \item{NoPg}{Number of progeny. A positive value creates that many lines; \code{-1} uses the allocation of an existing destination population.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of creating doubled-haploid progeny.
}

\details{
For each progeny, a parental individual is chosen randomly from the active individuals of the parental population and a haploid meiotic product is generated. The selected haploid chromosome set is then doubled so that both homologues of each progeny carry the same recombined haplotype.

The simulation genome must be diploid. For positive \code{NoPg}, that many doubled-haploid progeny are generated. The special value \code{-1} uses the allocation of an existing progeny population and therefore requires an existing destination with positive allocation.

The destination population is constructed from raw chromosome segments. Cached genotype and evaluation data can subsequently be recalculated. If meiosis or chromosome construction fails, the incomplete progeny population is removed.
}
