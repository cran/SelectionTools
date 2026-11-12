\name{gd.distance.similarity}
\alias{gd.distance.similarity}
\title{Gd distance similarity}

\description{
Calculates genetic distance or similarity measures, with optional bootstrap inputs.
}

\usage{
gd.distance.similarity(dta, measure, par, boot=0, no.pop.u=0, pop.u=0,
                       no.ind.u=0, ind.u=0, no.mar.u=0, mar.u=0)
}

\arguments{
   \item{dta}{Input data object.}
   \item{measure}{Distance or similarity measure.}
   \item{par}{Parameters used by the selected method.}
   \item{boot}{Bootstrap control.}
   \item{no.pop.u}{Number of populations used for bootstrap sampling.}
   \item{pop.u}{Population bootstrap indices.}
   \item{no.ind.u}{Numbers of individuals used for bootstrap sampling.}
   \item{ind.u}{Individual bootstrap indices.}
   \item{no.mar.u}{Number of markers used for bootstrap sampling.}
   \item{mar.u}{Marker bootstrap indices.}
}
\value{
A \code{data.frame} with columns \code{OTU1}, \code{OTU2}, and \code{Measure}, containing pairwise population comparisons. When resampling-based uncertainty is requested an additional \code{SDev} column contains the corresponding standard deviations.
}
