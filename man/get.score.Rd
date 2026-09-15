% file get.score.Rd
% copyright (C) 2001 Frisch/Maurer
%
\name{get.score}
\alias{get.score}
\title{get.score()}
\description{
  The function \emph{get.score()} returns a field with the calculated
  genotypic index for every individual of the population \emph{pop}. If
  an effect is specified then the value of this effect is returned.
}
\usage{
get.score(pop, effectfile=NULL)
}
\arguments{
   \item{pop}{Name of the evaluated population.}
   \item{effectfile}{Name of effect which should be evaluated.}
}
\value{
A one-column \code{data.frame} named \code{gvalue}, containing the genetic values obtained after genotyping and evaluating the requested population.
}
