% file genotype.population.Rd
% copyright (C) 2001 Frisch/Maurer
%
\name{genotype.population}
\alias{genotype.population}
\title{genotype.population()}
\description{
Calculates marker genotypes for one or more simulation populations from their chromosome-segment representation and the currently loaded linkage map.
}
\usage{
genotype.population(PopNames)
}
\arguments{
\item{PopNames}{One population name, a space-separated character string of population names, or a character vector of population names to genotype.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of generating/storing genotype information for the requested populations.
}

\details{
Simulation populations store chromosome segments as the primary genetic representation. This function projects those segments onto the loci of the currently loaded map and stores the two alleles for every active individual and mapped locus.

All active individuals of a population are treated as one genotype category: after a successful call every active individual has genotype data. If a population contains a mixture of genotyped and non-genotyped active individuals, the complete population is recalculated. Allocation failure clears the population genotype cache rather than leaving a partial mixture.

Population names may be supplied as a vector; the R wrapper combines multiple names into the space-separated form accepted by the C routine. Genotyping does not evaluate genetic values. Existing evaluations that depend on a genotype may need to be recalculated after changes to the underlying population or map.
}
