\name{st.load.vcf.data}
\alias{st.load.vcf.data}
\title{Loads Compact STvcf Marker Data with an Optional Linkage Map}
\description{
Loads a compact STvcf object directly into a SelectionTools marker data set
and, when present, its linkage map without intermediate marker or map text files.
}
\usage{
st.load.vcf.data(STvcf, data.set="default", strictACTG=FALSE)
}
\arguments{
\item{STvcf}{A list in STvcf format.}
\item{data.set}{Character string naming the SelectionTools marker data set to
be created or replaced.}
\item{strictACTG}{Logical, default \code{FALSE}. When enabled, the reader accepts only
A, C, G and T labels or the abbreviations Ade, Cyt, Gua and Thy, matched
case-insensitively, and assign fixed internal codes A/Ade=1, C/Cyt=2,
G/Gua=3 and T/Thy=4 independently of REF/ALT order. Other labels become
missing.}
}
\details{
With \code{strictACTG=TRUE}, nucleotide labels determine the internal codes.
Accepted labels are A, C, G, T, Ade, Cyt, Gua and Thy, case-insensitively.
Markers having none of these labels in REF or ALT are discarded.
A marker with an accepted label remains even if all its genotypes are
missing. A diploid call containing any excluded allele becomes entirely
missing. The surviving marker data and map remain aligned, and chromosome
codes are consecutive after filtering. Marker-only input remains marker-only.
If no marker remains, loading returns a failure status without replacing an
existing destination. The input STvcf object is not modified.

Strict mode ignores the optional \code{allele.codes} and
\code{prediction.state} attributes because labels define a new coding;
individual pedigree metadata is retained. For ordinary compact exchange of
SelectionTools data, including objects with numeric REF/ALT labels, retain
the default \code{strictACTG=FALSE}. It preserves the existing coding and
fitted-state behavior described below.

Fixed nucleotide coding supports matching markers with different REF/ALT
orders; it does not perform marker matching, strand complementation or
indel normalization.

The optional \code{allele.codes} attribute is a list of integer coding vectors, named and ordered by marker ID. It records the internal allele codes corresponding to REF and successive ALT entries. STvcf exchange retains this attribute so fixed or subsetted loci keep the same coding as their fitted marker reference. The loader validates it before replacing a destination. Conversions preserve it; a manually edited marker list must keep it aligned. It is not part of plain-text VCF exchange.

The optional versioned \code{pedigree} attribute is a data frame with \code{name} and \code{pedigree} columns and integer attribute \code{version=1L}. Its names must match the individual identifiers in order. Without it, individuals are founders. Full strings are owned independently of ancestor populations. The optional \code{prediction.state} attribute preserves a compatible active symbolic-A fitted state during STvcf/RDS exchange. Conversion in either direction retains both attributes; loading validates them before replacing the destination.

The optional \code{physical.map} attribute is a data frame with
\code{marker}, \code{chromosome} and \code{position} columns and
\code{position.unit = "bp"}. Its marker identifiers must be unique and
cover every imported marker; chromosome labels must be nonmissing and
positions finite and nonnegative. The loader copies this map into compiled
marker-data storage. It is independent of the genetic map in
\code{CHROM}/\code{POS}; a physical map never supplies cM distances.
Existing genetic-map inputs retain their previous interpretation.
\code{st.gwas} uses the stored physical coordinates when \code{map = NULL}.
Copying and subsetting marker data preserve the physical metadata, and
STvcf export includes it again. Export of physical-map-only data retains
missing genetic \code{CHROM}/\code{POS} and empty \code{CHROM.levels}.

STvcf is a compact exchange representation for diploid marker data with
VCF-style allele indexing and optional genetic positions. Loading
converts it into the existing SelectionTools marker and linkage-map data
structures; no parallel internal analysis representation is introduced.

The list must contain exactly the following nine components, each exactly
once. Extra, missing, duplicated, unnamed, or unknown list components are
rejected by the C loader:
\describe{
\item{\code{CHROM}}{Integer running chromosome codes, consecutive from 1, or all \code{NA_integer_} for marker-only data.}
\item{\code{POS}}{Double-precision genetic positions in cM. Values must be
finite and nonnegative for mapped data, or all \code{NA_real_} for marker-only data.}
\item{\code{ID}}{Unique non-empty marker identifiers containing fewer than
256 bytes. Whitespace and control characters are not permitted.}
\item{\code{REF}}{One reference-allele label per marker. REF may not be
missing, empty, \code{"."}, or comma-separated.}
\item{\code{ALT}}{\code{"."} or a comma-separated alternative-allele list.
ALT entries must be non-empty, unique, and different from REF.}
\item{\code{individual}}{Unique non-empty individual identifiers, subject to the R string-size limit. Whitespace and control characters are not permitted.}
\item{\code{allele1}}{Marker-by-individual raw or integer matrix containing
the first VCF allele index.}
\item{\code{allele2}}{Marker-by-individual raw or integer matrix containing
the second VCF allele index.}
\item{\code{CHROM.levels}}{Character labels associated with the running
chromosome codes. An empty character vector declares marker-only data.}
}

Individual names must be unique and non-empty, subject to the R string-size limit.
Whitespace and control characters are rejected; ordinary punctuation such as
dots, underscores, and hyphens is permitted.

Alleles use VCF-style indexing: \code{0} denotes the reference allele
(\code{REF}), and \code{1}, \code{2}, and subsequent indices refer to
the alleles listed in \code{ALT}. When loaded into SelectionTools,
these are represented by allele codes \code{1}, \code{2}, \code{3},
and so on.

The input object can be prepared with \code{\link{st.dataframe.to.STvcf}}.
This function converts a data frame of marker information and
genotypes into the STvcf format required here.

Each marker must have one reference allele. Alternative alleles must
be unique and differ from the reference allele. \code{"."} denotes
\code{ALT} when no alternative allele is defined; empty entries and
dots within a list of alternative alleles are not allowed.

Markers may have up to 254 alleles, excluding missing values.
Multiallelic markers can be loaded, but genomic prediction methods
requiring biallelic markers must be used with a suitable subset.
Fixed loci are accepted by the symbolic-A relationship path and contribute zero; other methods retain their own marker requirements.

Only diploid genotypes are supported. If either allele is missing,
the entire genotype is treated as missing. The order of non-missing
alleles is preserved. The logical \code{phased} attribute indicates
whether the supplied genotypes are phased. If this attribute is absent,
the data are treated as unphased. The attribute must be logical and
must not contain missing values. Loading does not infer phase.

Mapped input must provide genetic map positions in centimorgans for all markers.
Markers are ordered by chromosome and position. Markers sharing the
same position retain their relative input order, but their positions
are adjusted slightly to make them distinct. The supplied STvcf object
is not modified. \code{\link{st.get.map}} reports the map
positions used by SelectionTools.

For marker-only input without a linkage map, every \code{CHROM} entry is \code{NA_integer_}, every \code{POS} entry is \code{NA_real_}, and \code{CHROM.levels} is \code{character()}. All nine components remain required. Partial maps, nonmissing positions without chromosome levels, and incomplete chromosome assignments are rejected. Marker-only loading preserves input marker order and installs no map. It supports marker statistics and relationship-matrix construction; map-dependent operations need a genuine linkage map. Further details appear in \code{\link{v-flintmaize-vcf}} for supplied marker-only data.

Marker statistics are calculated from the supplied data. An existing
data set named \code{data.set} is replaced only if loading succeeds.
Otherwise, it remains unchanged. Performance data are not loaded by
this function.
}
\value{
An integer status, invisibly: zero on success and a negative value on rejected input. The supplied STvcf marker and map information is installed in \code{data.set} only after successful staging.
}

\examples{
data("v-tropmaize-vcf", package = "SelectionTools")
marker.frame <- st.STvcf.to.dataframe(v.tropmaize.vcf)[1:12, c(1:6, 7:14)]
vcf <- st.dataframe.to.STvcf(marker.frame)
st.load.vcf.data(vcf, data.set = "maize.example")
st.get.map(data.set = "maize.example")
}

\seealso{
\code{\link{st.dataframe.to.STvcf}},
\code{\link{st.markerdata.to.STvcf}},
\code{\link{st.STvcf.to.dataframe}}.
}
