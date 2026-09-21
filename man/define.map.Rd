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
   \item{description}{A string containing a description of the linkage map}
}
\value{
An invisible \code{list} returned by \code{generate.map.file()} or \code{linkage.map.load()}. Its \code{retval} component is the status code of the map-generation/loading step; the main effect is defining the simulation map.
}

\details{
The compact description consists of comma-separated locus specifications with
the grammar
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
currently active simulation linkage-map definition only after the description
has been accepted by the map-generation/loading steps.
}
