\name{remove.all.populations}
\alias{remove.all.populations}
\title{remove.all.populations()}

\description{
Removes all currently supported simulation populations.
}

\usage{
remove.all.populations()
}
\value{
Returns invisible NULL. The main effect is removal of all populations.
}

\details{
The current population list is collected and all listed populations are removed. Genome parameters, linkage-map definitions, effect definitions, and marker-data sets are separate objects and are not removed with the simulation populations.
}

\seealso{
For vector based input of loci, see \code{\link{evaluate.allele.freq2}}.
}