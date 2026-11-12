\name{st.detect.core}
\alias{st.detect.core}
\title{Detect a Haplotype Core Collection}
\description{
Selects compact core collections from one phased marker data set using either
forward selection or exact best-subset selection on distinct block-specific
haplotypes.
}
\usage{
st.detect.core(pool, n.max,
               haplotype.options = NULL, blocks = NULL,
               exclude.compatible = TRUE, method = "forward")
}
\arguments{
\item{pool}{Name of the phased marker data set containing the candidate genotypes.}
\item{n.max}{Positive integer maximum number of genotypes in the core.}
\item{haplotype.options}{Named list of options passed to \code{st.def.hblocks}:
\code{hap}, \code{hap.unit}, \code{ld.threshold}, \code{ld.criterion},
\code{tolerance}, and \code{hap.symbol}. An additional \code{ld.measure}
option selects the measure for \code{st.calc.ld} when LD-based blocks are
requested. Exactly one of \code{haplotype.options} and \code{blocks} is required.
An empty list requests the default five-marker blocks.}
\item{blocks}{Block-definition data frame, typically returned by
\code{st.def.hblocks}. The \code{Name}, \code{Chrom}, and \code{Markers}
columns define unique block names, chromosomes, and semicolon-separated marker
membership. Blocks must be nonoverlapping and contain consecutive markers on
one chromosome.}
\item{exclude.compatible}{Logical. \code{TRUE} selects the conservative
missing-compatible rule. In forward selection, a complete haplotype is not
counted as a marginal gain when it is compatible with an incomplete haplotype
already present in the selected core. In best-subset selection, a represented
complete haplotype is excluded from the optimization score when it is compatible
with an incomplete haplotype anywhere in that candidate subset.}
\item{method}{Selection method. \code{"forward"} adds one genotype at a time
using the largest current marginal gain. \code{"best.subset"} finds the exact
optimal subset separately for each core size up to \code{n.max}.}
}
\details{
The marker data set must already have consistent gametic phase. Phasing is not
performed by this function. Haplotype identity is the ordered sequence of
marker allele states within a block. The same sequence in two different blocks
counts as two block-specific haplotypes. Only complete haplotype sequences are
counted as represented segments.

With \code{method="forward"}, selection starts with an empty core. The first
genotype is the one representing the largest number of distinct complete
block-specific haplotypes. At every subsequent step, the remaining genotype
contributing the largest number of previously unrepresented qualifying
haplotypes is added. Shared haplotypes count once. Equal gains are resolved by
bytewise lexicographic genotype name.

With \code{method="best.subset"}, an exact branch-and-bound search is carried
out separately for core sizes 1 through \code{n.max}. For each size, the subset
maximizing the number of qualifying represented haplotypes is returned. When
\code{exclude.compatible = TRUE}, a complete variant is excluded from this
objective if it is compatible with an incomplete haplotype in the candidate
subset. Branch-and-bound bounds ignore future compatibility masking only when
forming an upper bound; because additional masking can only lower a descendant
score, this pruning remains exact. Equal optimal subsets are resolved by the
first combination of genotype names in bytewise lexicographic order. Runtime
can grow combinatorially for large pools or large values of \code{n.max}; no
heuristic result is silently substituted.

With \code{exclude.compatible = TRUE}, incomplete haplotypes in the already
selected core can make otherwise new complete variants compatible with the
current core. Such variants are excluded from the marginal selection gain.
The \code{segments} values nevertheless report the cumulative number of
distinct complete block-specific haplotypes actually represented by the
selected genotypes.

Forward selection stops when \code{n.max} genotypes have been selected or when
no remaining genotype has a positive marginal gain. Exact best-subset search
stops early once an optimal subset represents every complete block-specific
haplotype observed in the pool. Thus either method can return fewer than
\code{n.max} core sizes.

Counting and both search methods are implemented in C. OpenMP builds use the
existing \code{st.set.num.threads} setting. Exact search parallelizes independent
first-genotype branches; forward search parallelizes candidate scoring. Results,
including ties, are independent of thread count. Input marker data and block
definitions are not modified, and no random numbers are consumed.
}
\value{
A list with two elements:
\describe{
\item{\code{genotypes}}{A list of core sets. With forward selection, element
\code{k} contains all genotype names selected through step \code{k}. With
best-subset selection, element \code{k} is the exact optimal subset of size
\code{k}; its names are in bytewise lexicographic order.}
\item{\code{segments}}{An integer vector parallel to \code{genotypes} giving
the cumulative number of distinct complete block-specific haplotypes represented
by each core set.}
}
If no complete haplotype can be represented, both elements have length zero.
}
\seealso{\code{st.select.enrichment}, \code{st.calculate.enrichment},
\code{st.def.hblocks}, \code{st.phase}}
