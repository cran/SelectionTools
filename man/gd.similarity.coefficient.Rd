\name{gd.similarity.coefficient}

\alias{gd.similarity.coefficient}

\title{gd.similarity.coefficient()}

\description{
Calculates matching coefficients between individuals.
}

\usage{
gd.similarity.coefficient(dta, measure="jac", par="dist", boot=0, no.pop.u=0,
                          pop.u=0, no.ind.u=0, ind.u=0, no.mar.u=0, mar.u=0,
                          resample.primers=FALSE)
}
\arguments{
   \item{dta}{Name of the dataset containing the observed marker
   genotypes. The format of the dataset is illustrated by the dataset
   gd.salad.aflp.}
   \item{measure}{\code{"dic"} is the Dice coefficient,
   \code{"jac"} the Jaccard similarity, and  
   \code{"sma"} is the simple matching coefficient.
   Default is the Jaccard similarity.}
   \item{par}{Parameters to be calculated:
   \tabular{ll}{
     \code{"dist"} \tab Genetic similarity measure  \cr
     \code{"sdev/jackm"} \tab Jackkniving over markers \cr
     \code{"sdev/bootm"} \tab Bootstrapping over markers \cr
     }}
   \item{boot}{Number of bootstrap runs for estimating the standard
   deviation of the  estimates.}
   \item{no.pop.u}{Function argument.}
   \item{pop.u}{Function argument.}
   \item{no.ind.u}{Function argument.}
   \item{ind.u}{Function argument.}
   \item{no.mar.u}{Function argument.}
   \item{mar.u}{Function argument.}
   \item{resample.primers}{if set to \code{TRUE} Resampling is done over
   primers, instead of markers. Usefull for AFLP data. Primers must be
   separated in the dataset from the markers by a \code{+}.}
}

\note{
Each population should consist of one indivdiual only.
  }

\value{
A \code{data.frame} with columns \code{OTU1}, \code{OTU2}, and \code{Measure}, containing the requested pairwise similarity coefficients; an additional \code{SDev} column is included when resampling-based standard deviations are requested. Returns \code{NULL} invisibly on validation failure.
}
