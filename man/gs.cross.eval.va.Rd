\name{gs.cross.eval.va}
\alias{gs.cross.eval.va}
\title{Evaluates Crosses by Segregation Variance}

\description{
Calculates the expected segregation variance for every possible pair of
parents in a SelectionTools data set from previously estimated marker
effects.
}

\usage{
gs.cross.eval.va(data.set = "default", pop.type = "unlinked", t = 0,
                 map.function = "Haldane")
}

\arguments{
  \item{data.set}{Character string naming the SelectionTools data set. The
  data set must contain estimated marker effects.}
  \item{pop.type}{Character string specifying the progeny population type.
  The implementation accepts \code{"unlinked"}, \code{"DH"}, \code{"SSD"},
  and \code{"CRS"}.}
  \item{t}{For \code{"DH"} and \code{"SSD"}, the nonnegative integer number of
  additional random-intermating generations before line development. Zero
  corresponds to F1-derived DH lines or F2-derived SSD lines. For SSD,
  additional intermating starts from the F2; this is not the number of
  selfing generations. The parameter is method-specific for other linked
  population types.}
  \item{map.function}{Character string specifying the map function used to
  obtain recombination frequencies when linkage is taken into account. The
  default is \code{"Haldane"}; \code{"Kosambi"} is also supported.}
}

\details{
The function estimates the segregation variance \eqn{\sigma_g^2}{sigma_g^2}
of each cross from the marker effects. The segregation variance describes how
strongly progeny from a cross are expected to differ genetically.

For \code{pop.type = "DH"} or \code{"SSD"}, linkage and recombination are
taken into account using the genetic map; this is the cross-variance
calculation described by Osthushenrich et al. (2017). With
\code{pop.type = "unlinked"}, the calculation treats the loci as unlinked.
The current implementation also accepts \code{"CRS"} as a linked population
type; its existing calculation uses one generation irrespective of \code{t}.


For two inbred parents, \code{pop.type = "SSD"} describes the limiting
variance among completely homozygous lines obtained by repeated selfing.
For two loci with recombination frequency \eqn{r}, the parental-origin
correlation in these lines is
\deqn{\rho_{SSD}(r,t)=\frac{1-2r}{1+2r}(1-r)^t.}
The corresponding probability of matching parental origin is
\eqn{q=(1+\rho_{SSD})/2}. At \code{t = 0}, this gives
\eqn{q=1/(1+2r)}. Different chromosomes are independent. Finite selfing
generations are represented by \code{ssd.mating}; they are not specified
by increasing \code{t} in the analytical variance calculation.

The arguments are scalar values. Invalid population types, map functions or
nonfinite, negative or fractional values of \code{t} produce an error before
the stored cross scores change.

The calculated variance is stored as the cross criterion \code{va}. Its
square root, \eqn{\sigma_g}{sigma_g}, is used by \code{gs.cross.eval.es} when
the expected superior progeny criterion is calculated. Results for all
evaluated crosses can subsequently be obtained with \code{gs.cross.info}.
}

\references{
Osthushenrich T, Frisch M, Herzog E (2017) Genomic selection of crossing
partners on basis of the expected mean and variance of their derived lines.
PLoS ONE 12(12):e0188839.
}

\value{
On a successful run, an invisible implementation-level \code{list} returned
by the final timing/information helper; the meaningful result is the side
effect of evaluating all crosses by segregation variance.
}

\seealso{
\code{gs.cross.eval.mu}, \code{gs.cross.eval.es}, \code{gs.cross.info}
}
