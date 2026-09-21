\name{gs.cross.eval.gd}
\alias{gs.cross.eval.gd}
\title{Evaluate Crosses by Genetic Distance}

\description{
Calculates a genetic distance for every possible pair of individuals in a
SelectionTools data set and stores it as a cross-evaluation criterion.
}

\usage{
gs.cross.eval.gd(dist = "rd", data.set = "default")
}

\arguments{
  \item{dist}{Character string specifying the genetic-distance measure.
  SelectionTools uses \code{"rd"} for Rogers distance, \code{"mrd"} for
  Modified Rogers distance, and \code{"euc"} for Euclidean distance.}
  \item{data.set}{Character string naming the SelectionTools data set.}
}

\details{
Genetic distance can be evaluated independently of estimated marker effects.
The result is stored as the cross criterion \code{gd} and can subsequently be
obtained together with the other cross criteria using \code{gs.cross.info}.

The definitions and formulas for Euclidean distance, Rogers distance, and
Modified Rogers distance are documented in \code{st.genetic.distances}. The
properties of these and related similarity and dissimilarity coefficients in
plant breeding are discussed by Reif et al. (2005).
}

\references{
Reif JC, Melchinger AE, Frisch M (2005) Genetical and mathematical properties
of similarity and dissimilarity coefficients applied in plant breeding and
seed bank management. Crop Science 45:1--7.
}

\value{
On a successful run, an invisible implementation-level \code{list} returned
by the final timing/information helper; the meaningful result is the side
effect of evaluating all crosses by genetic distance. Validation failures may
return \code{NULL}.
}

\seealso{
\code{st.genetic.distances}, \code{gs.cross.eval.mu}, \code{gs.cross.info}
}
