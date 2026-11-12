\name{st.calculate.enrichment}
\alias{st.calculate.enrichment}
\alias{st.select.enrichment}
\title{Calculates Haplotype Enrichment and Selects Donor Genotypes}
\description{
Counts additional block-specific haplotype alleles supplied by one phased
marker data set to another, or selects donor genotypes by exact best-subset
search or stepwise forward selection. Both stored homologues are considered.
}
\usage{
st.calculate.enrichment(base.pool, enrichment.pool,
                        haplotype.options = NULL, blocks = NULL,
                        exclude.compatible = TRUE)

st.select.enrichment(base.pool, enrichment.pool, n,
                     haplotype.options = NULL, blocks = NULL,
                     exclude.compatible = TRUE, method = "best.subset")
}
\arguments{
\item{base.pool}{Name of the already phased base marker data set.}
\item{enrichment.pool}{Name of the already phased donor marker data set.
For calculation, all its genotypes contribute. For selection, it is the
candidate pool.}
\item{n}{Positive integer number of donor genotypes to select.}
\item{haplotype.options}{Named list of options passed to
\code{st.def.hblocks}: \code{hap}, \code{hap.unit}, \code{ld.threshold},
\code{ld.criterion}, \code{tolerance}, and \code{hap.symbol}.
An additional \code{ld.measure} option selects the measure for
\code{st.calc.ld} when LD-based blocks are requested (default \code{"r2"}).
Omitted entries use the defaults of those routines. An empty list requests
the default five-marker blocks. Exactly one of
\code{haplotype.options} and \code{blocks} is required.}
\item{blocks}{Block-definition data frame, typically returned by
\code{st.def.hblocks}. The \code{Name}, \code{Chrom} and \code{Markers}
columns define unique block names, chromosomes, and semicolon-separated
marker membership. Blocks must be nonoverlapping and each must contain
consecutive markers on one chromosome. A subset of blocks is allowed.
Additional columns in the supplied table are retained.}
\item{method}{Selection method: \code{"best.subset"} (default) finds an
exact optimal subset; \code{"forward"} adds the genotype with the largest
remaining gain at each step. Forward selection is faster for large pools but
need not find the globally best subset.}
\item{exclude.compatible}{Logical. If \code{TRUE}, exclude additional
variants compatible with an incomplete base haplotype from the selection
objective and reported \code{Score}. Both categories are always reported
by \code{st.calculate.enrichment}.}
}
\details{
The two data sets must already have consistent gametic phase along each
homologue. Phasing is not performed by these functions. The data sets must
contain the same marker identities, chromosome assignments and map positions.
Genotypes are aligned by marker name. Their internal integer allele coding
must have the same biological meaning at each marker in both data sets.
SelectionTools does not retain the original textual REF/ALT labels in marker
data, so matching biological coding cannot be verified automatically. Creating
both pools by copying and restricting a common marker data set preserves the
coding. Individual names are local to each pool and may overlap between pools.

When options are supplied, block definitions are computed exclusively from a
private working copy of the base population, including its LD if needed.
The same base data and options therefore produce the same blocks when the
full donor pool is replaced by selected donors. Supplying fixed blocks instead
uses their marker membership directly. The input marker data sets, including
their block definitions, are not modified. The temporary compiled working data
set is released on return or error; no new simulation or phasing code is involved.

Haplotype identity is the ordered sequence of marker allele states within a
block, not an independently assigned haplotype code. Only complete sequences
are counted as observed alleles. A variant absent from complete base calls
is potentially hidden if it agrees with every observed marker of at least
one incomplete base haplotype. An entirely missing base haplotype is compatible
with every variant at that block. This is conditional on the supplied phase
and calls; lack of a compatible incomplete sequence does not prove biological
novelty. compiled allele arrays are read directly. Missing values on one
homologue do not discard observed values on the other homologue.

Selection maximizes the union of qualifying additional alleles, not the sum
of individual donor scores. Shared variants count once per block. With
\code{method="best.subset"}, an exact
branch-and-bound search prunes only branches that cannot improve the current
best score. Equal optima are resolved by the first combination of genotype
names sorted in bytewise lexicographic order. The result is deterministic and
does not consume random numbers. Runtime can grow combinatorially for large
pools; no heuristic result is silently substituted for an exact optimum.

With \code{method="forward"}, each step chooses the remaining genotype
contributing the most qualifying alleles not already supplied by selected
genotypes. Equal gains are resolved by bytewise lexicographic genotype name.
Selection continues to exactly \code{n} genotypes even when remaining gains
are zero. Block borders and compatibility with incomplete base haplotypes
remain fixed relative to the original base population throughout the search.
This makes both methods optimize the same enrichment objective.

Counting and both search methods are implemented in C. OpenMP builds use the
existing \code{st.set.num.threads} setting for parallel search; selections,
including ties, are independent of thread count. Serial builds use the same
algorithms. No random numbers are consumed.
}
\value{
\code{st.select.enrichment} returns a character vector of exactly \code{n}
donor genotype names: sorted by name for best-subset search, and in selection
order for forward search. It does not create a selected marker data set.

\code{st.calculate.enrichment} returns a list containing:
\describe{
\item{\code{summary}}{One-row data frame with \code{Base.lines},
\code{Enrichment.lines}, \code{Main}, \code{Additional}, \code{Hidden},
\code{No.match}, \code{Score}, \code{Enriched}, and \code{Gaining.blocks}.
\code{Additional = Hidden + No.match}; \code{Enriched = Main + Additional}.
Counts are sums of distinct block-specific alleles, not allele copies.}
\item{\code{by.block}}{The same allele counts and score for each block,
including blocks with zero gain.}
\item{\code{alleles}}{Additional observed variants, with columns
\code{Block}, \code{Sequence} and logical \code{Hidden}.}
\item{\code{blocks}}{The supplied or generated block-definition table.}
\item{\code{exclude.compatible}}{The scoring option used.}
}
}
\seealso{\code{st.phase}, \code{st.def.hblocks}, \code{st.calc.ld},
\code{st.copy.marker.data}, \code{st.restrict.marker.data}}
