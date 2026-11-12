\name{set.eff.weight}
\alias{set.eff.weight}
\title{set.eff.weight()}

\description{
  Assigns a weight to a defined effect
}

\usage{
set.eff.weight(fname, weight)
}

\arguments{
   \item{fname}{Name of the effect}
   \item{weight}{Weight of the effect}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of setting the requested effect weight.
}

\note{
   Effect weights are used when \code{evaluate.genotype} is
   carried out \emph{without} specifing one effect. In this case all loaded
   effects are used to determin tne genotypic value, and weights are
   applied.
   If evaluate genotype is used \emph{with} specifiying an effect, the
   weights are \emph{not} used. 
   \cr
   Further details appear in \code{define.effects} for an example

 }
