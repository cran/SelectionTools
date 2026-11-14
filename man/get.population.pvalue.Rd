\name{get.population.pvalue}
\alias{get.population.pvalue}
\title{Returns phenotypic values from a simulation population}

\description{
Returns stored phenotypic values for the active individuals of a simulation
population.
}

\usage{
get.population.pvalue(name, EffName=NULL)
}

\arguments{
\item{name}{Character string naming the population.}
\item{EffName}{Optional character string naming an effect.}
}

\details{
Phenotypic values must previously have been generated with
\code{phenotype.population}. If \code{EffName} is omitted, element 0 of the
phenotypic value array is returned. This is the weighted phenotypic selection
index. If \code{EffName} is supplied, the corresponding effect-specific
phenotypic value is returned.
}

\value{
A one-column \code{data.frame} named \code{pvalue}, containing one phenotypic value per active individual in the requested population.
}

\seealso{
\code{phenotype.population}, \code{get.population.gvalue}
}
