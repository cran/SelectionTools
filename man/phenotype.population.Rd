\name{phenotype.population}
\alias{phenotype.population}
\title{Generate phenotypic values for a simulation population}

\description{
Generates phenotypic values for all loaded effects and all active individuals
of a simulation population.
}

\usage{
phenotype.population(PopName, hsq)
}

\arguments{
\item{PopName}{Character string naming one simulation population.}
\item{hsq}{Numeric vector containing one heritability for each loaded effect,
in the order returned by \code{list.effects()}. Each value must be greater
than zero and not greater than one.}
}

\details{
Phenotypic values are generated for all loaded effects in one operation.
Partial phenotypic evaluation is not supported. The length of \code{hsq}
must therefore equal the number of loaded effects.

The function first genotypes and evaluates the population. For effect
\eqn{i}, let \eqn{V_{G,i}} be the sample variance of the effect-specific
genetic values and let \eqn{h_i^2} be the supplied heritability. The residual
variance is
\deqn{V_{E,i} = V_{G,i}/h_i^2 - V_{G,i}.}
For every individual an independent normal residual with variance
\eqn{V_{E,i}} is added to the effect-specific genetic value.

The phenotypic value array is parallel to the genetic value array. Elements
1 through the number of effects contain the effect-specific phenotypic
values. Element 0 contains the weighted phenotypic selection index,
calculated with the current effect weights in the same way as the genetic
selection index.

Phenotypic values are population-wide derived data. Re-evaluating genetic
values, explicitly setting genetic values, removing an evaluation, removing
cached genotypes, changing the effect configuration or effect weights, or
removing the linkage map invalidates stored phenotypic values.
}

\value{
On successful phenotype generation, an invisible implementation-level \code{list} returned by the compiled routine; its \code{retval} component is the status code. Validation failures return \code{NULL} invisibly.
}

\seealso{
\code{evaluate.population}, \code{get.population.pvalue},
\code{population.sort}, \code{select.n.best}, \code{select.all.best}
}
