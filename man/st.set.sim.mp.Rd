\name{st.set.sim.mp}
\alias{st.set.sim.mp}
\title{Sets sim mp used by SelectionTools}

\description{
Copies the mapped loci of a SelectionTools marker data set into the active simulation linkage map.
}

\usage{
st.set.sim.mp(data.set="default")
}

\arguments{
\item{data.set}{Character string naming the SelectionTools marker data set that provides the map points.}
}
\value{
An invisible \code{list} returned by the compiled routine. Component \code{r} is the status value used by the simulation bridge; the function is primarily called for its side effect of transferring the requested data to the simulation backend.
}

\details{
Map points are transferred with chromosome, position, marker name, and class information. The transfer respects the SelectionTools distinction between physical map storage and map access order and preserves the source access relationship used by the marker matrix.

Simulation genome parameters should already be compatible with the marker-data map. This stage installs map points but does not create the simulation population.
}
