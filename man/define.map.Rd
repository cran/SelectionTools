\name{define.map}
\alias{define.map}
\title{define.map()}

\description{
Defines and installs a simulation linkage map from a compact textual map description.
}

\usage{
define.map(description)
}

\arguments{
   \item{description}{A character string or character vector containing
   locus specifications. Multiple specifications can be separated by
   commas or whitespace within a string, or supplied as separate
   elements of a character vector.}
}

\value{
Invisibly returns a list whose \code{retval} component indicates
success or failure. The function is primarily called for its side
effect of defining the simulation linkage map.
}

\details{
Genome parameters must be defined before calling this function, for
example using \code{\link{set.genome.par}}. Each call replaces the
currently defined linkage map.

\preformatted{[*Number] Locusname Classname Position}
where \code{*Number} is optional. \code{Locusname} is the locus-name prefix,
\code{Classname} is the locus class, and \code{Position} controls placement.

For a single explicitly positioned locus, \code{Position} is written as
\code{chromosome/position}. The chromosome is one-based and the position must
lie between 0 and the chromosome length defined for the simulation genome.
The position is in the same unit as the simulation chromosome lengths
(Morgan in the standard SelectionTools simulation setup). A multiplicity
prefix is not allowed for explicitly positioned loci.

Generated positions can be specified as follows:
\describe{
  \item{\code{random}}{Generate \code{Number} loci at independently random
  positions over the complete genome. Chromosomes are sampled in proportion
  to their lengths, so these positions are genuinely random rather than
  equidistant.}
  \item{\code{terminal}}{Distribute \code{Number} loci equally among the
  chromosomes and place them equidistantly on each chromosome including the
  telomeres. \code{Number} must be a multiple of the number of chromosomes.}
  \item{\code{nonterminal}}{Distribute \code{Number} loci equally among the
  chromosomes and place them equidistantly in the chromosome interiors,
  excluding the telomeres. \code{Number} must be a multiple of the number of
  chromosomes.}
}
For generated positions a running integer suffix is added to the supplied
\code{Locusname} prefix to obtain unique locus names.

The description is first converted to a temporary map file by
\code{generate.map.file} and is then loaded through the simulation
linkage-map loader. The generated map is validated for chromosome bounds,
position bounds, unique locus names, and complete specifications, and is
ordered by chromosome and position before loading. The operation replaces the
currently active simulation linkage map definition only after the description
has been accepted by the map-generation/loading steps.
}

\examples{
set.genome.par(no.chrom = 1,
               no.hom = 2,
               chrom.len = 0.05)

# Define a single locus
# Format: Name Class Chromosome/Position
define.map("form trait1 1/0")

# Replace the map with two loci supplied in one string
define.map("form trait1 1/0, color trait2 1/0.02")

# Define the same map using a character vector
define.map(c("form trait1 1/0",
             "color trait2 1/0.02"))

# Inspect the resulting map
linkage.map.get()
}