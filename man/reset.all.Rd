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
Returns invisble \code{NULL}. The main effect is resetting all package simulation state handled by the function.
}
\note{
All informations are lost, besides the genome parameter
(\code{genome.parameter.set}).
}
\details{
Storage owned by simulation populations is returned to the package's compiled
memory pool. The pool retains backing blocks for reuse by subsequent
simulations. Further details appear in \code{population.remove} for the storage lifecycle.
}
