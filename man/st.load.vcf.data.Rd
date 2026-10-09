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

Individual names must be unique, contain only letters and digits (no dots, hyphens, etc.),
and have between 1 and 255 characters.

Alleles use VCF-style indexing: \code{0} denotes the reference allele
(\code{REF}), and \code{1}, \code{2}, and subsequent indices refer to
the alleles listed in \code{ALT}. When loaded into SelectionTools,
these are represented by allele codes \code{1}, \code{2}, \code{3},
and so on.

Prepare the input object using \code{\link{st.dataframe.to.STvcf}}.
This function converts a data frame of marker information and
genotypes into the STvcf format required here.

Each marker must have one reference allele. Alternative alleles must
be unique and differ from the reference allele. Use \code{"."} for
\code{ALT} when no alternative allele is defined; empty entries and
dots within a list of alternative alleles are not allowed.

Markers may have up to 254 alleles, excluding missing values.
Multiallelic markers can be loaded, but genomic prediction methods
requiring biallelic markers must be used with a suitable subset.
Monomorphic markers must also be excluded before using these methods.

Only diploid genotypes are supported. If either allele is missing,
the entire genotype is treated as missing. The order of non-missing
alleles is preserved. The logical \code{phased} attribute indicates
whether the supplied genotypes are phased. If this attribute is absent,
the data are treated as unphased. The attribute must be logical and
must not contain missing values. Loading does not infer phase.

All markers must have genetic map positions in centimorgans.
Markers are ordered by chromosome and position. Markers sharing the
same position retain their relative input order, but their positions
are adjusted slightly to make them distinct. The supplied STvcf object
is not modified. Use \code{\link{st.get.map}} to inspect the map
positions used by SelectionTools.

Marker statistics are calculated from the supplied data. An existing
data set named \code{data.set} is replaced only if loading succeeds.
Otherwise, it remains unchanged. Performance data are not loaded by
this function.
}
\value{
\code{NULL}, returned invisibly. The function is called for its side effect of loading the supplied STvcf marker and map information into \code{data.set}.
}

\examples{
# Three markers and three individuals.
# POS contains genetic positions in cM.
# 0 denotes REF and 1 denotes ALT.
markers <- data.frame(
   CHROM = c(1, 1, 2),
   POS = c(5, 15, 10),
   ID = c("M1", "M2", "M3"),
   REF = c("A", "C", "G"),
   ALT = c("G", "T", "A"),
   Plant1 = c("0/0", "0/1", "1/1"),
   Plant2 = c("0/1", "1/1", "0/0"),
   Plant3 = c("1/1", "0/0", "0/1")
)

# Convert the data frame to the compact STvcf format
marker.data <- st.dataframe.to.STvcf(markers)

# Load marker genotypes and their map into SelectionTools
st.load.vcf.data(marker.data, data.set = "example")

# Inspect the loaded linkage map
mp <- st.get.map(data.set = "example")
mp
}

\seealso{
\code{\link{st.dataframe.to.STvcf}},
\code{\link{st.markerdata.to.STvcf}},
\code{\link{st.STvcf.to.dataframe}}.
}