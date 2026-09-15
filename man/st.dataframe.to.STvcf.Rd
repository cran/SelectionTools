\name{st.dataframe.to.STvcf}
\alias{st.dataframe.to.STvcf}
\title{Convert a VCF-Like Data Frame to Compact STvcf Format}
\description{
Converts a wide diploid VCF-like R data frame to the compact STvcf list
representation used by \code{st.load.vcf.data()}.
}
\usage{
st.dataframe.to.STvcf(data)
}
\arguments{
\item{data}{A data frame containing VCF-style fixed fields and sample
genotypes. The required layout and conversion conventions are described
below.}
}
\details{
The converter is deliberately strict. It is intended to create only STvcf
objects that can subsequently be accepted by the SelectionTools marker-data
and simulation-transfer path. Invalid marker or individual names, malformed
genotypes, inconsistent REF/ALT definitions, invalid map positions, and
unsupported allele indices are rejected during conversion rather than being
left for a later consumer to discover.

The first five data-frame columns are fixed by position and must be, in this
order, \code{CHROM} or \code{#CHROM}, \code{POS}, \code{ID}, \code{REF},
and \code{ALT}. This positional convention is intentional. Sample names are
allowed to be ordinary alphanumeric strings such as \code{POS}, \code{INFO},
or \code{FORMAT}; therefore fixed fields cannot safely be identified only by
searching all column names.

Data-frame column names must be non-missing and non-empty. Duplicate names are
permitted only when they arise because a sample has the same name as a fixed
or metadata field. Sample names themselves must be unique.

Marker IDs are mandatory and unique. Marker and individual identifiers may
contain only ASCII letters \code{A-Z} and \code{a-z} and digits
\code{0-9}. Dots, underscores, hyphens, whitespace, punctuation, and
non-ASCII characters are rejected. Names must contain fewer than 256 bytes.
The same restrictions are checked again by \code{st.load.vcf.data()} and by
the reverse marker-data conversion. Names are never silently modified.

Input chromosome labels may be numeric or character. They are converted to
consecutive running integer codes 1, 2, ..., in order of first appearance.
The original input labels are retained in \code{CHROM.levels}. One
information message at level 0 reports the mapping. Chromosome labels are not
marker or individual names and are not subject to the alphanumeric name
restriction.

STvcf positions are genetic positions in centimorgan. If a valid metadata
column named \code{CM} is supplied, it is used for the genetic position and
the input \code{POS} values are not used for mapping. Otherwise
\code{POS} is assumed to contain cM. This is important for conventional VCF
data, where POS normally denotes physical base-pair position.

When a valid \code{FORMAT} column is present, columns between ALT and FORMAT
may be the VCF metadata fields \code{QUAL}, \code{FILTER}, \code{INFO}, and
the SelectionTools-specific \code{CM} field. These metadata fields must be
unique. Every column after FORMAT is a sample column, irrespective of its
name. This allows sample identifiers such as \code{POS}, \code{INFO},
\code{CM}, or \code{FORMAT} without ambiguity.

A FORMAT column is recognized only when every marker-specific definition
contains exactly one non-empty GT field. Empty FORMAT fields and duplicated GT
entries are invalid. The FORMAT definition may vary among markers. The GT
subfield is extracted from every sample cell. Other subfields, including DP,
AD, GQ, and phase-set fields, are ignored and are not retained. A non-missing
sample field that is too short to contain the GT subfield defined by FORMAT,
or that contains an empty GT subfield, is rejected as malformed rather than
being converted to missing. If FORMAT is absent, columns after the first five
fixed fields are direct GT sample columns, except that an alphanumeric column
named CM whose complete contents
are valid finite nonnegative numeric values is interpreted as the genetic CM
field. This rule is the convention used to distinguish an optional CM field
from a sample that itself is named CM.

Genetic positions are stored as R double-precision numeric values without
intentional rounding. Factor CM or POS columns are converted through their
printed character values before numeric conversion, so factor level numbers
are never used as map positions. Positions must be finite and nonnegative.

REF must contain exactly one non-missing allele label. It may not be
\code{"."} and may not contain a comma. ALT may be \code{"."}, denoting no
defined alternative allele, or a comma-separated list. Every ALT entry must
be non-empty, may not be \code{"."}, must differ from REF, and must be unique
within the marker. Thus definitions such as REF=A with ALT=A or ALT=C,C are
rejected.

Multiallelic markers are accepted. VCF allele index 0 denotes REF, index 1
the first ALT allele, index 2 the second ALT allele, and so on. Every observed
non-missing allele index must be present in the corresponding REF/ALT
definition. At most 253 ALT alleles may be defined for a marker. Together
with REF this corresponds to at most 254 biological allele states in a
SelectionTools marker data set; missing data use a separate internal state.

Only diploid genotypes are supported. Both \code{"/"} and \code{"|"} are
accepted. The supplied allele order is retained in \code{allele1} and
\code{allele2}, but the separator itself is not stored. Retaining allele order
is not a statistical phasing procedure and must not be interpreted as phase
inference. Phase-set boundaries and other phase metadata are not represented
in STvcf.

Accepted complete-missing representations include R \code{NA}, an empty
field, \code{"."}, \code{"./."}, \code{".|."}, \code{"-1"},
\code{"-1/-1"}, \code{"-1|-1"}, and the character value \code{"NA"}.
The \code{-1} spelling is accepted as a SelectionTools backward-compatibility
extension; standard VCF output continues to use \code{"."} for missing
alleles. Partially missing diploid genotypes such as \code{0/.}, \code{./1},
\code{0|.}, \code{.|1}, \code{0/-1}, and \code{-1/1} are normalized to
complete missingness. The non-missing member of a partial pair must
nevertheless be a syntactically valid nonnegative integer allele index within
the REF/ALT definition.
Malformed fields such as an empty allele component or a nonnumeric
non-missing component are rejected rather than being hidden by the
partial-missing rule.

The converter writes \code{allele1} and \code{allele2} as raw matrices.
Valid non-missing VCF allele indices are 0 through 253, and raw value 255
denotes missing data. Integer allele matrices are also accepted by the STvcf
loader and reverse converter for manually constructed STvcf objects, but the
253-ALT contract still applies.

The returned list has class \code{"STvcf"} and components \code{CHROM},
\code{POS}, \code{ID}, \code{REF}, \code{ALT}, \code{individual},
\code{allele1}, \code{allele2}, and \code{CHROM.levels}. Marker and
individual order are retained. Marker ordering is changed only when the STvcf
object is installed in the internal SelectionTools linkage-map structure.

The conversion is not a byte-for-byte VCF round trip. QUAL, FILTER, INFO,
FORMAT subfields other than GT, phase-set information, physical POS when CM
is supplied, and the distinction between slash and vertical-bar separators
are not stored.
}
\value{
A \code{list} of class \code{"STvcf"} containing marker metadata, individual identifiers, two compact allele matrices, and the chromosome-label mapping. Positions are stored as double-precision values.
}
