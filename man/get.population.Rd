\name{get.population}
\alias{get.population}
\title{get.population()}

\description{
Returns the raw chromosome-segment representation of a simulation population as an R data frame.
}

\usage{
get.population(PopName)
}

\arguments{
\item{PopName}{Character string naming the simulation population.}
}

\details{
Each row identifies an individual, chromosome, homologue, segment start position, and allele. These rows are the primary simulation representation from which marker genotypes are derived when \code{genotype.population} is run.

The returned data frame can be used to inspect or serialize a population and is the representation accepted by \code{init.population}. Retrieving the population does not alter genotype or evaluation caches.
}
\value{
A \code{data.frame} with columns \code{ind}, \code{chrom}, \code{hom}, \code{pos}, and \code{all}, representing the active population in the package population format.
}
