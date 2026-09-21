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
The simulation routines and the marker-data/genomic-prediction routines use
different internal representations. Marker data sets store genotypes in
marker matrices, whereas the simulation backend reconstructs individual
chromosomes and uses a count-location representation of chromosome segments
for the simulation of meiosis. This conversion is therefore required before
a SelectionTools marker data set can be used by the simulation routines.

The function first calls \code{reset.all}. Consequently, previously defined
simulation populations, effects, and map information are discarded before the
new simulation state is initialized. It then performs the bridge stages in
sequence: \code{st.set.sim.gp}, \code{st.set.sim.mp},
\code{st.set.sim.ef}, and \code{st.set.sim.pp}. Processing stops immediately
if a stage reports failure.

The source data set must provide marker data and a valid linkage map for the
genome-parameter and map stages. The chromosome parameters, linkage map, and
marker genotypes are transferred to the simulation backend and used to
reconstruct the simulation chromosomes.

Estimated marker effects are optional. If no estimated marker effects are
present, \code{st.set.sim.ef} succeeds without installing an effect and
population transfer continues normally. If effects are present, they are
transferred directly under \code{effect.name}; no intermediate effect file is
created. The genomic-selection routines normally use additive marker dosage
coding 0/1/2, whereas the simulation representation is equivalent to -1/0/1.
The required baseline correction of the general mean is therefore applied
during effect transfer. See \code{st.set.sim.ef} for the exact coding and
formula.

The individual bridge functions can also be called separately. In particular,
a simulation may be initialized once with genome parameters, map, and optional
effects and then populated from additional compatible marker data sets using
repeated calls to \code{st.set.sim.pp}. This staged route avoids the
\code{reset.all} performed by another call to \code{st.set.simpop}.
}

\value{
Invisibly returns the \code{list} produced by the last completed simulation-bridge stage. Component \code{r} is the bridge status value. The main effect is initialization of the simulation backend from the SelectionTools data set.
}

\references{
Maurer HP, Melchinger AE, Frisch M (2008) Population genetic simulation and
data analysis with Plabsoft. Euphytica 161:133--139.

Stranden I, Christensen OF (2011) Allele coding in genomic evaluation.
Genetics Selection Evolution 43:25. See Equation (2).
}

\seealso{
\code{st.set.sim.gp}, \code{st.set.sim.mp}, \code{st.set.sim.ef},
\code{st.set.sim.pp}, \code{reset.all}
}
