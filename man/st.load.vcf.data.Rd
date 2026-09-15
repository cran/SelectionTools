\name{st.load.vcf.data}
\alias{st.load.vcf.data}
\title{Load Compact STvcf Marker and Linkage-Map Data}
\description{
Loads a compact STvcf object directly into a SelectionTools marker data set
and installs its linkage map without intermediate marker or map text files.
}
\usage{
st.load.vcf.data(STvcf, data.set="default")
}
\arguments{
\item{STvcf}{A list in STvcf format.}
\item{data.set}{Character string naming the SelectionTools marker data set to
be created or replaced.}
}
\details{
STvcf is a compact exchange representation for diploid marker data with
VCF-style allele indexing and one genetic position per marker. Loading
converts it into the existing SelectionTools marker and linkage-map data
structures; no parallel internal analysis representation is introduced.

The list must contain exactly the following nine components, each exactly
once. Extra, missing, duplicated, unnamed, or unknown list components are
rejected by the C loader:
\describe{
\item{\code{CHROM}}{Integer running chromosome codes, consecutive from 1.}
\item{\code{POS}}{Double-precision genetic positions in cM. Values must be
finite and nonnegative.}
\item{\code{ID}}{Unique marker identifiers containing only ASCII letters and
digits and fewer than 256 bytes.}
\item{\code{REF}}{One reference-allele label per marker. REF may not be
missing, empty, \code{"."}, or comma-separated.}
\item{\code{ALT}}{\code{"."} or a comma-separated alternative-allele list.
ALT entries must be non-empty, unique, and different from REF.}
\item{\code{individual}}{Unique individual identifiers containing only ASCII
letters and digits and fewer than 256 bytes.}
\item{\code{allele1}}{Marker-by-individual raw or integer matrix containing
the first VCF allele index.}
\item{\code{allele2}}{Marker-by-individual raw or integer matrix containing
the second VCF allele index.}
\item{\code{CHROM.levels}}{Character labels associated with the running
chromosome codes.}
}

Name validation is repeated in C. A manually assembled STvcf therefore cannot
bypass the rule that marker and individual names contain only letters and
digits. Dots, underscores, hyphens, whitespace, punctuation, non-ASCII
characters, empty names, duplicate names, and names of 256 bytes or more are
rejected.

STvcf allele index 0 denotes REF and indices 1, 2, ... denote ALT alleles. On
loading these become SelectionTools biological allele codes 1, 2, 3, ... .
At most 253 ALT alleles may be defined at one marker, so the largest resulting
SelectionTools biological allele code is 254. The two allele matrices must use
the same storage type. Raw STvcf matrices use 255 for missing; integer matrices
use \code{NA_integer_} for missing data. Missingness is converted to the
SelectionTools internal missing representation.

REF and ALT definitions are validated independently of the R converter. A
REF value repeated in ALT, a duplicated ALT value, an empty ALT token, or a
dot embedded inside an ALT list is rejected. Thus malformed manually created
STvcf objects cannot create distinct internal integer alleles for identical
textual allele labels.

The SelectionTools marker representation allows at most 254 biological
allele codes per marker, plus the internal missing state. Multiallelic markers
are valid input as long as they fit this representation. The ordinary statistical genomic
selection functions that construct one design-matrix column per marker are
biallelic; multiallelic data must be restricted to an appropriate biallelic
marker set before those functions are called. Monomorphic markers are also
unsuitable for that ordinary biallelic design-matrix construction.

Only diploid data are supported. If either stored allele is missing, the
complete internal diploid genotype is set to missing. The order of two non-missing alleles is
retained, and those two arrays become the two homologues when marker data are
transferred to simulation. This does not imply statistical phase inference.

All STvcf markers require map positions. Markers are loaded in chromosome and
cM order. Exact position ties retain STvcf row order and are separated by
small deterministic internal offsets because the existing linkage-map and
simulation structures require distinct ordered positions. The STvcf object
itself is not changed; \code{st.get.map()} reports the internally used map
positions.

Normal SelectionTools marker statistics are constructed after conversion.
They are calculated from exactly the markers present in the STvcf object.
Compact raw or integer R storage is expanded into the existing C genotype and
marker-statistics structures during loading.

The operation replaces \code{data.set} transactionally. The replacement is
committed only after names, alleles, genotype matrices, linkage-map data, and
marker statistics have been validated and constructed successfully. On
failure, the previous contents of the named marker data set remain intact.
Performance data are not installed by this function.
}
\value{
\code{NULL}, returned invisibly. The function is called for its side effect of loading the supplied STvcf marker and map information into \code{data.set}.
}
