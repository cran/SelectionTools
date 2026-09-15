\name{remove.all.populations}
\alias{remove.all.populations}
\title{remove.all.populations()}

\description{
Removes all currently registered simulation populations.
}

\usage{
remove.all.populations()
}
\value{
Returns the same invisible implementation-level \code{list} as \code{population.remove.all()}; the main effect is removal of all populations.
}

\details{
The current population list is collected and all listed populations are removed. Genome parameters, linkage-map definitions, effect definitions, and marker-data sets are separate objects and are not removed with the simulation populations.
}

