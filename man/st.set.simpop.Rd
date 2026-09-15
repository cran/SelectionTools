\name{st.set.simpop}
\alias{st.set.simpop}
\title{Initialize the Simulation Backend from a SelectionTools Data Set}

\description{
Initializes the simulation backend from a SelectionTools marker data set and
creates a diploid simulation population.
}

\usage{
st.set.simpop(pop.name, data.set="default", effect.name="effect")
}

\arguments{
\item{pop.name}{Character string naming the simulation population to create.
The name must contain 2 to 255 ASCII letters or digits.}
\item{data.set}{Character string naming the source SelectionTools marker data
set.}
\item{effect.name}{Character string naming transferred marker effects in the
simulation backend. The default is \code{"effect"}.}
}

\details{
The function first calls \code{reset.all} and then performs the bridge stages
in sequence: \code{st.set.sim.gp}, \code{st.set.sim.mp},
\code{st.set.sim.ef}, and \code{st.set.sim.pp}. Processing stops immediately
if a stage reports failure.

The source data set must provide marker data and a valid linkage map for the
genome-parameter and map stages. Estimated marker effects are optional. If no
estimated marker effects are present, \code{st.set.sim.ef} succeeds without
installing an effect and population transfer continues normally.

If effects are present, the intercept and marker effects are transferred
directly to the simulation backend under \code{effect.name}; no intermediate
effect file is created.

The individual bridge functions can also be called separately. In particular,
a simulation may be initialized once with genome parameters, map, and optional
effects and then populated from additional compatible marker data sets using
repeated calls to \code{st.set.sim.pp}.
}

\value{
Invisibly returns the \code{list} produced by the last completed simulation-bridge stage. Component \code{r} is the bridge status value. The main effect is initialization of the simulation backend from the SelectionTools data set.
}

\seealso{
\code{st.set.sim.gp}, \code{st.set.sim.mp}, \code{st.set.sim.ef},
\code{st.set.sim.pp}
}
