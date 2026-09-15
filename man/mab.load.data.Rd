\name{mab.load.data}
\alias{mab.load.data}

\title{Loads the data of a simulation}

\description{
  Loads a simulation result from the explicitly supplied \code{result.file}
  for further analysis or printing. Relative file names are resolved against
  \code{st.data.dir}; absolute paths are used directly.}

\usage{
mab.load.data(simulation.name, simulation.run, result.file)
}

\arguments{
   \item{simulation.name}{Simulation-name identifier retained for compatibility with the MAB workflow. It does not determine the file name.}
   \item{simulation.run}{Simulation-run identifier retained for compatibility with the MAB workflow. It does not determine the file name.}
   \item{result.file}{Required name of the result file. Relative names are resolved against \code{st.data.dir}; absolute paths are used directly. No file name is generated from \code{simulation.name} or \code{simulation.run}.}
}
\value{
The object \code{outp} loaded from the requested simulation result file. For files created by \code{mab.simulate()}, this is the named simulation-result \code{list} described there.
}

\seealso{
  \code{\link{mab.simulate}}, \code{\link{mab.tabulate}},
  \code{\link{mab.Examples}}
}
