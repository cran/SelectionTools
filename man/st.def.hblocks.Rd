% copyright (C) 2001 Frisch/Maurer
%
\name{st.def.hblocks}
\alias{st.def.hblocks}
\title{Defines Haplotype Blocks}
\description{
Defines haplotype-block boundaries by a fixed number of markers, a genetic
window in cM, linkage disequilibrium, or disjunct chromosome segments.
}
\usage{
st.def.hblocks(hap = 5, hap.unit = 1, ld.threshold = 0.8,
               ld.criterion = "none", tolerance = 0, hap.symbol = "b",
               out.filename = "blocks", auxfiles = TRUE, data.set = "default")
}
\arguments{
   \item{hap}{Block size for the fixed-block methods: number of markers when
   \code{hap.unit=1}, or genetic distance in cM when \code{hap.unit=2}.}
   \item{hap.unit}{Block-definition mode when \code{ld.criterion="none"}:
   \code{1} builds blocks by marker count, \code{2} builds blocks over a
   specified genetic distance in cM, and \code{0} builds disjunct chromosome
   segments as used for introgression populations.}
   \item{ld.threshold}{Minimum LD value used when \code{ld.criterion} is
   \code{"average"} or \code{"flanking"}. Previously calculated LD values
   must be available in \code{data.set}.}
   \item{ld.criterion}{Criterion used for LD-based block definition.
   \code{"none"} uses \code{hap.unit}; \code{"average"} evaluates the
   average pairwise LD within a candidate block, and \code{"flanking"}
   evaluates LD between the markers flanking a candidate block.}
   \item{tolerance}{Non-negative number of neighboring markers used as a
   tolerance region when extending an LD-defined block. Within that region
   the best value of the selected LD criterion is used for the extension
   decision.}
   \item{hap.symbol}{Character prefix used for generated block names. With the
   default \code{"b"}, blocks are named \code{b000000}, \code{b000001}, and
   so on.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Character string naming the marker data set in which
   haplotype blocks are defined.}
}

\details{
Haplotype analysis is performed in two steps. This function defines the block
boundaries along each chromosome; \code{st.recode.hil} subsequently searches
the marker sequences between those boundaries and recodes the observed
haplotype variants as alleles of the block loci.

For \code{hap.unit=1}, consecutive observed markers are grouped into blocks of
up to \code{hap} markers without crossing chromosome boundaries. For
\code{hap.unit=2}, consecutive markers are grouped according to their genetic
map positions in cM. For \code{hap.unit=0}, block boundaries are determined
from disjunct marker segments in the current genotypes.

If \code{ld.criterion} is \code{"average"} or \code{"flanking"}, the
LD-based method takes precedence over \code{hap.unit}. It requires LD values
to have been calculated beforehand. The algorithm starts blocks from adjacent
marker pairs exceeding \code{ld.threshold} and extends them while the selected
LD criterion remains at least the threshold. Markers that cannot be joined to
another marker under the criterion form single-marker blocks.

This function only defines the blocks; it does not replace the marker data
by haplotype loci. That replacement occurs when \code{st.recode.hil} is
called.
}

\value{
When \code{auxfiles = TRUE} and the operation succeeds, invisibly returns a
\code{data.frame} describing the defined haplotype blocks:
\describe{
  \item{\code{Chrom}}{Chromosome number of the block.}
  \item{\code{Pos}}{Block position, calculated by the current implementation
  as the midpoint between the positions of the first and last marker in the
  block.}
  \item{\code{Name}}{Generated haplotype-block name.}
  \item{\code{Class}}{Class assigned to the haplotype block by the compiled
  implementation.}
  \item{\code{Markers}}{Names of the markers used to construct the block,
  stored as a semicolon-separated marker list in the intermediate
  representation.}
}
Returns \code{NULL} when no auxiliary result is requested or when the
operation does not produce a readable block description.
}

\seealso{
\code{st.calc.ld}, \code{st.recode.hil}
}
%
