\name{population.info.get}
\alias{population.info.get}
\alias{population.pedigree}
\title{Returns stored information for an individual}

\description{
Returns the character information stored for an individual in a simulation population.
}

\usage{
population.info.get(PopName, ind)

population.pedigree(name, source = c("auto", "simulation", "marker"))
}

\arguments{
   \item{PopName}{Name of the population.}
   \item{ind}{Individual index or indices.}


\item{name}{One simulation-population or marker-dataset name.}
\item{source}{Registry to query. The default requires an unambiguous existing name.}

}
\value{
A character value containing the stored information for the requested individual when the individual index is valid; otherwise \code{NULL}.
}

\section{Whole-population simulated pedigree}{
Returns stable individual names and their complete stored simulated pedigree expressions from a simulation population or marker dataset.

Imported founders retain their exact imported identifiers. New crosses combine the stored pedigree expressions of the actual parents as \code{(left_x_right)}; DH and SSD production append the respective operation inside parentheses. The complete expression is stored when the individual is created and remains available after ancestor populations are removed. Copying, moving, sorting, selection and population renaming retain each individual's identity. This accessor does not infer a pair of immediate parents from a complex pedigree.

The result follows current active individual order. Repeated copies in a simulation population may have the same identity. Export to a marker dataset rejects duplicated individual identifiers. A missing query does not create a marker dataset. A name present in both registries requires an explicit \code{source}.

A fresh data frame with character columns \code{name} and \code{pedigree}, using automatic row names. Empty valid populations return zero rows.

}
