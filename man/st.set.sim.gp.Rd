\name{st.set.sim.gp}
\alias{st.set.sim.gp}
\title{Sets sim gp used by SelectionTools}

\description{
Derives simulation genome parameters from a SelectionTools marker-data map and installs them in the simulation backend.
}

\usage{
st.set.sim.gp(data.set="default")
}

\arguments{
\item{data.set}{Character string naming the SelectionTools marker data set that provides the linkage map.}
}
\value{
An invisible \code{list} returned by the compiled routine. Component \code{r} is the status value used by the simulation bridge; the function is primarily called for its side effect of transferring the requested data to the simulation backend.
}

\details{
The marker data set must contain a valid linkage map. Chromosome identifiers must define a contiguous sequence from chromosome 1 through the maximum chromosome number, and positions must be finite and non-negative. Chromosome lengths are derived from the mapped positions.

The SelectionTools marker matrix is diploid, so the bridge installs genome parameters compatible with two homologues. This is the genome-parameter stage of transferring a marker data set into the simulation representation. It does not by itself install map points or create a population.
}
