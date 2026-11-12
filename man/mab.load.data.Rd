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
The object \code{outp} loaded from the requested simulation result file. For
files created by \code{mab.simulate()}, this is a named list containing the
simulation parameters and, depending on the requested analyses, the following
summaries:
\describe{
  \item{\code{TargetAlleles}}{Target-allele summary for selected plants over
  generations.}
  \item{\code{RPG.gw}}{Recurrent-parent genome content over the whole genome.}
  \item{\code{RPG.reg1}, \code{RPG.reg2}, \dots}{Recurrent-parent genome
  content for separately defined chromosome regions.}
  \item{\code{DonorSegments.gw}}{Donor chromosome segments in selected
  individuals, counted per homologue.}
  \item{\code{DonorSegments.reg1}, \code{DonorSegments.reg2}, \dots}{Donor
  chromosome segments in separately defined chromosome regions.}
  \item{\code{LinkageDrag1}, \code{LinkageDrag2}, \dots}{Donor chromosome
  segment attached to each target locus; target loci are evaluated separately
  and linkage drag is averaged over the two homologues.}
  \item{\code{SM}}{Number of marker data points attributed to single-marker
  assays.}
  \item{\code{HT}}{Number of marker data points attributed to high-throughput
  assays.}
}
Further details appear in \code{mab.simulate} for the selection-strategy rules used to count
\code{SM} and \code{HT} marker data points.
}

\seealso{
  \code{\link{mab.simulate}}, \code{\link{mab.tabulate}},
  \code{\link{mab.Examples}}
}
