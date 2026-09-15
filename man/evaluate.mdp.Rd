\name{evaluate.mdp}
\alias{evaluate.mdp}
\title{Evaluate marker-assisted donor-parent information for populations}

\description{
Evaluate marker-assisted donor-parent information for populations.
}

\usage{
evaluate.mdp(NamePg, NameP1, NameP2, effectfile=NULL)
}

\arguments{
   \item{NamePg}{Name of the population to evaluate.}
   \item{NameP1}{Name of the first parental population.}
   \item{NameP2}{Name of the second parental population.}
   \item{effectfile}{Effect file or effect specification.}
}
\value{
When a specific \code{effectfile} is supplied, the result is the invisible implementation-level \code{list} returned by the compiled MDP evaluator. When \code{effectfile = NULL}, the function evaluates all active effects for side effects and returns \code{NULL}.
}
