\name{st.get.simpop.perfdata}
\alias{st.get.simpop.perfdata}
\title{Return simulation-population data in SelectionTools format}

\description{
Imports a simulation population as marker data and uses an evaluated simulation effect as the phenotype vector of the SelectionTools data set.
}

\usage{
st.get.simpop.perfdata(pop.name, data.set, effect)
}

\arguments{
\item{pop.name}{Character string naming the simulation population.}
   \item{data.set}{Character string naming the target SelectionTools marker data set.}
   \item{effect}{Effect name or effect specification passed to population evaluation.}
}
\value{
An invisible integer status value returned by \code{gs.set.performance.data()} after simulation-population performance data have been transferred to the SelectionTools data set.
}

\details{
The function first imports the simulation population into \code{data.set}. It then ensures that the simulation population is genotyped, evaluates the requested effect, retrieves the resulting population genetic values, and writes those values into the phenotype vector of the imported SelectionTools marker data set.

The operation therefore combines a simulation-to-marker-data transfer with generation of performance data from the simulation effect model. The diploidy and naming requirements of the simulation-population importer apply. In particular, the population name used by this transfer must contain only ASCII letters and digits, because \code{st.get.simpop()} constructs alphanumeric marker-data individual identifiers from that name.
}
