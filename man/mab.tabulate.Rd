\name{mab.tabulate}
\alias{mab.tabulate}

\title{Tabulate the results of a series of simulations}

\description{
  Tabulates results from a series of simulations. The files to load are
  supplied explicitly in \code{result.files}.}

\usage{
mab.tabulate(simulation.names.runs, par.summary = TRUE, ext.summary = FALSE,
             RPG.gw.mean = FALSE, TargetAlleles = FALSE, RPG.gw.Q10 = FALSE,
             RPG.reg.Q10 = FALSE, RPG.reg.mean = FALSE, DonorSegments.gw = FALSE,
             DonorSegments.reg = FALSE, LinkageDrag = FALSE, MarkerAnalyses = FALSE,
             digits = 1, percent = TRUE, result.files)
}

 \arguments{
   \item{simulation.names.runs}{Character vector containing alternating
    simulation names and simulation-run identifiers. It is interpreted as
    a two-column matrix, one row per simulation result to be tabulated.}
   \item{par.summary}{Logical. Include a summary of the main simulation parameters.}
   \item{ext.summary}{Logical. Include an extended summary of the simulation parameters.}
   \item{RPG.gw.mean}{Logical. Include mean genome-wide recurrent-parent genome (RPG) values.}
   \item{TargetAlleles}{Logical. Include the mean target-allele summary.}
   \item{RPG.gw.Q10}{Logical. Include the Q10 value of the genome-wide RPG.}
   \item{RPG.reg.Q10}{Logical. Include Q10 values of RPG in the defined chromosome regions.}
   \item{RPG.reg.mean}{Logical. Include mean RPG values in the defined chromosome regions.}
   \item{DonorSegments.gw}{Logical. Include the mean number of donor segments genome-wide.}
   \item{DonorSegments.reg}{Logical. Include the mean number of donor segments in the defined chromosome regions.}
   \item{LinkageDrag}{Logical. Include mean linkage-drag summaries for the target loci.}
   \item{MarkerAnalyses}{Logical. Include mean single-marker (SM) and high-throughput (HT) marker-analysis counts.}
   \item{digits}{Number of decimal places used when rounding numerical
    summaries.}
   \item{percent}{Logical. If \code{TRUE}, RPG summaries are multiplied by
    100; if \code{FALSE}, they remain on the proportion scale.}
   \item{result.files}{Character vector containing one explicit result-file name for each scenario in \code{simulation.names.runs}, in the same order. Relative names are resolved against \code{st.data.dir}; absolute paths are used directly.}
}
\value{
A named \code{list} of the summary matrices requested by the logical arguments. Possible components include \code{par.summary}, \code{ext.summary}, \code{TargetAlleles.mean}, \code{RPG.gw.mean}, \code{RPG.gw.Q10}, numbered regional \code{RPG.reg*.mean} and \code{RPG.reg*.Q10} matrices, \code{DonorSegments.gw.mean}, numbered \code{DonorSegments.reg*.mean} matrices, numbered \code{LinkageDrag*.mean} matrices, \code{SM.mean}, and \code{HT.mean}. If incompatible scenario definitions are detected, \code{NULL} is returned invisibly.
}

\seealso{
  \code{\link{mab.simulate}}, \code{\link{mab.load.data}},
  \code{\link{mab.compare}}, \code{\link{mab.Examples}}
}
