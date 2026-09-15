% file get.population.gvalue.Rd
% copyright (C) 2001 Frisch/Maurer
%
\name{get.population.gvalue}
\alias{get.population.gvalue}
\title{Get genetic values from a simulation population}

\description{
Returns stored evaluated genetic values for the active individuals of a
simulation population.
}

\usage{
get.population.gvalue(name, EffName=NULL)
}

\arguments{
\item{name}{Character string naming the population.}
\item{EffName}{Optional character string identifying the effect-specific
values to retrieve.}
}

\details{
Genetic values are derived population data created by population evaluation
or explicitly supplied through \code{set.population.gvalue}.

If \code{EffName} is omitted, element 0 of the genetic value array is
returned. After ordinary population evaluation this is the weighted genetic
selection index. If \code{EffName} is supplied, the corresponding
effect-specific genetic value is returned.
}

\value{
A one-column \code{data.frame} named \code{gvalue}, containing one genetic value per active individual in the requested population.
}

\seealso{
\code{evaluate.population}, \code{set.population.gvalue},
\code{get.population.pvalue}
}
