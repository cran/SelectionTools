\name{st.STvcf.to.dataframe}
\alias{st.STvcf.to.dataframe}
\title{Converts Compact STvcf Data to a VCF-Like Data Frame}
\description{
Expands a compact STvcf list to a wide VCF-like R data frame with one marker
per row and one genotype column per individual.
}
\usage{
st.STvcf.to.dataframe(STvcf)
}
\arguments{
\item{STvcf}{A compact STvcf list containing marker metadata, individual
identifiers, chromosome mapping, and two allele matrices.}
}
\details{
The optional \code{allele.codes} attribute is a list of integer coding vectors, named and ordered by marker ID. It records the internal allele codes corresponding to REF and successive ALT entries. STvcf exchange retains this attribute so fixed or subsetted loci keep the same coding as their fitted marker reference. The loader validates it before replacing a destination. Conversions preserve it; a manually edited marker list must keep it aligned. It is not part of plain-text VCF exchange.

The optional versioned \code{pedigree} attribute is a data frame with \code{name} and \code{pedigree} columns and integer attribute \code{version=1L}. Its names must match the individual identifiers in order. Without it, individuals are founders. Full strings are owned independently of ancestor populations. The optional \code{prediction.state} attribute preserves a compatible active symbolic-A fitted state during STvcf/RDS exchange. Conversion in either direction retains both attributes; loading validates them before replacing the destination.

The function validates the STvcf object before expansion. The list must
contain exactly the nine documented STvcf components, each exactly once;
extra, missing, duplicated, unnamed, or unknown components are rejected.
Component types and dimensions must agree, chromosome codes must be
consecutive and used, map
positions must be finite and nonnegative, allele indices must agree with
REF/ALT, and marker and individual identifiers must obey the same restrictions
as \code{st.load.vcf.data()}.

Marker identifiers must be unique, non-empty, and shorter than 256 bytes. Individual identifiers must be unique and non-empty; they use dynamically allocated strings subject to the R string-size limit. Whitespace and control characters are rejected; ordinary
punctuation such as dots, underscores, and hyphens is permitted. These checks
are intentionally the same on conversion in both directions.

REF may not be missing, \code{"."}, empty, or comma-separated. ALT may be
\code{"."} or a comma-separated list. ALT entries must be non-empty, unique,
and different from REF. Invalid or duplicated allele labels are rejected.

The returned data frame contains \code{CHROM}, \code{POS}, \code{ID},
\code{REF}, \code{ALT}, and \code{FORMAT}, followed by one sample column per
entry in \code{STvcf$individual}. FORMAT is written as \code{"GT"} for every
marker. Including FORMAT is deliberate: it makes the sample boundary
unambiguous even when an individual identifier is itself an alphanumeric name
such as POS, INFO, CM, or FORMAT.

CHROM is written using \code{CHROM.levels}. POS is written as the current
STvcf double-precision genetic position in cM without intentional rounding.

Genotypes use VCF allele indices, with 0 for REF and 1, 2, ... for ALT
alleles. Allele order is retained. If the STvcf object has logical attribute
\code{phased=TRUE}, non-missing genotypes are written with \code{"|"} and
missing genotypes as \code{".|."}. Otherwise \code{"/"} and \code{"./."}
are used. STvcf objects created before this attribute was introduced, or
manually created objects without it, are treated as unphased.

The two allele matrices must use the same storage type. Raw matrices use raw
value 255 for missing data; integer matrices use \code{NA_integer_} for
missing data. STvcf objects accepted by SelectionTools may define at most 253
ALT alleles per marker, so valid non-missing VCF allele indices range from 0
through 253. Allele indices are checked against the corresponding REF/ALT
definition.

STvcf does not retain QUAL, FILTER, INFO, non-GT FORMAT fields, phase-set
metadata, or a physical POS field that was superseded by CM. These fields
cannot be reconstructed. For an STvcf object produced by
\code{st.dataframe.to.STvcf()}, conversion to a data frame and back preserves
the data-set level phased/unphased state through the logical \code{phased}
attribute.
}
\value{
A VCF-like \code{data.frame} with fixed columns \code{CHROM}, \code{POS}, \code{ID}, \code{REF}, \code{ALT}, and \code{FORMAT}, followed by one genotype column for each individual.
}
