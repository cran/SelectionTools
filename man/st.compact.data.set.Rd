% copyright (C) 2001 Frisch/Maurer
%
\name{st.compact.data.set}
\alias{st.compact.data.set}
\alias{st.repack.data.set}
\title{Compacts a Marker Data Set and Repacks Allele Storage}
\description{
Releases calculated results and reusable workspace allocations associated with
one SelectionTools marker data set while retaining its marker data and linkage
map.
\code{st.repack.data.set()} selects compact allele storage while retaining
the data set and calculated results.
}
\usage{
st.compact.data.set(data.set = "default")

st.repack.data.set(data.set = "default", commit = TRUE)
}
\arguments{
  \item{data.set}{Character string naming the marker data set.}
  \item{commit}{Logical. For repacking, TRUE replaces storage; FALSE verifies and reports the proposed representation without changing the data set.}
}
\details{
The retained core consists only of individual and marker names, the allele storage containing the marker genotypes, the data-set phase state, and the
linkage map including its marker-to-data indices.

All other dynamically allocated data-set state is discarded. This includes
marker statistics, phenotypic and genotypic values, effects and model results,
design and relationship matrices, haplotype-block and cross structures,
chromosome-specific recombination/LD workspaces, and reusable temporary
buffers. Consequently columns such as \code{effects}, \code{phenotypic.value}, and \code{genotypic.value} in
\code{st.data.set.list()} can change from \code{TRUE} to \code{FALSE} after
compaction.

\code{st.repack.data.set()} retains calculated state and chooses the smallest
supported representation from the maximum number of observed nonmissing alleles
at any marker: at most three gives four bits per diploid genotype, four to
fifteen gives eight bits, and sixteen or more gives sixteen bits. Four- and
eight-bit storage includes a per-marker dictionary that preserves original
allele codes, missing values, and homologue order. Dictionary memory is additional
to the genotype payload. Imports and genotype replacements select this storage
automatically. Repacking is useful after changes that reduce allele diversity.
New storage is fully constructed and checked before replacing live storage.

Values that are required again must be recalculated or reloaded by the
corresponding SelectionTools function.
The dedicated \code{genomic_relationship_A} cache and dependent selected
hybrid SCA caches are cleared. Hybrid performance and parental references
remain supported, and requested relationship matrices are rebuilt on demand.
}
\value{
Invisibly returns the estimated amount of work memory released, in MiB. The
estimate excludes memory allocator overhead.
\code{st.repack.data.set()} returns a list with \code{bits.per.genotype},
\code{genotype.bytes}, \code{dictionary.bytes}, \code{roundtrip.verified},
and \code{committed}. The byte counts describe the proposed or committed
allele representation, excluding names, maps, and allocator overhead.
}
