% file reset.all.Rd
% copyright (C) 2001 Frisch/Maurer
%
\name{reset.all}
\alias{reset.all}
\title{reset.all()}
\description{
  "Resets" the program into the initial state.
}
\usage{
reset.all()
}
\value{
Returns the same invisible implementation-level \code{list} as \code{reset.mdp()} after resetting the other simulation state. The main effect is resetting all package simulation state handled by the function.
}
\note{
All informations are lost, besides the genome parameter
(\code{genome.parameter.set}).
}
