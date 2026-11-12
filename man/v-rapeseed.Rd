\name{v-rapeseed}
\alias{v-rapeseed}
\alias{v.rapeseed.vcf}
\alias{v.rapeseed.phe}
\alias{v.rapeseed.map}
\docType{data}
\title{Rapeseed Marker and Flowering-Time Data for the GWAS Vignette}
\description{
Rapeseed data for genome-wide association analysis, containing 22857 SNP markers
for 193 rapeseed lines and flowering-time observations, provided by
Rod Snowdon (University of Giessen). Genotypes use the compact STvcf representation;
physical marker coordinates and phenotypes are separate objects.
}
\usage{
data("v-rapeseed", package = "SelectionTools")
}
\format{
The data file loads three objects:
\describe{
\item{\code{v.rapeseed.vcf}}{A nine-component list of class
\code{"STvcf"} containing 22857 markers and 193 individuals. The
\code{ID} and \code{individual} components retain the supplied unique
marker and line identifiers in their original order. \code{REF} and
\code{ALT} contain the neutral allele labels \code{"1"} and \code{"2"}.
The \code{allele1} and \code{allele2} components are 22857 by 197 raw
matrices using VCF allele indices 0 and 1. All calls are complete.
\code{CHROM} is entirely \code{NA_integer_}, \code{POS} is entirely
\code{NA_real_}, and \code{CHROM.levels} is \code{character()}, declaring
marker-only data without a genetic linkage map. The \code{physical.map}
attribute contains the physical map described below, with positions in
base pairs and rows aligned exactly to \code{ID}. The \code{phased}
attribute is \code{FALSE}.}
\item{\code{v.rapeseed.phe}}{A data frame with 193 rows and columns
\code{gen} (genotype identifier) and \code{flt} (flowering time). Its row order matches \code{v.rapeseed.vcf$individual}.
All 193 flowering-time values are observed. Lines \code{SW068},
\code{SW125}, \code{SW146}, and \code{SW195}, which lack phenotypes,
are excluded from both marker and phenotype data. The measurement
unit and time origin are not specified.}
\item{\code{v.rapeseed.map}}{A data frame with 22857 rows and columns
\code{marker}, \code{chromosome}, and \code{position}. It preserves the
source marker identifiers, numeric chromosome labels 1 through 19, and
integer physical positions in base pairs. The \code{position.unit}
attribute is \code{"bp"}. Rows match the marker order in the STvcf
object.}
}
}
\details{
All 22857 markers and the 193 individuals with observed flowering time
are retained in the packaged data. The four lines without flowering-time
observations are excluded before packaging. No marker filtering,
imputation, phasing, or phenotype replacement is performed.

The original genotype scores \code{-1}, \code{0}, and \code{1} are
represented by unphased VCF calls \code{0/0}, \code{0/1}, and \code{1/1},
respectively. Numeric REF/ALT labels are used because nucleotide
identities are not available. With the default
\code{strictACTG = FALSE}, \code{st.load.vcf.data} maps these allele
indices to internal SelectionTools allele codes 1 and 2. The resulting
default marker dosage is the original score plus one. No sign or allele
frequency information is lost.

The source positions are physical coordinates, whereas STvcf map fields
are reserved for genetic positions in centimorgans. Consequently the
STvcf object contains no linkage map. Loading
the STvcf object with \code{st.load.vcf.data} copies its \code{physical.map}
attribute to compiled marker-data storage. Subsequent
\code{st.gwas} calls use these coordinates automatically, without a
separate \code{map} argument.
The separately loaded \code{v.rapeseed.map} is identical to this attribute.
The loader does not interpret the physical-map attribute as a genetic map.
These data do not supply a genetic map for recombination or simulation.

Two within-chromosome coordinate ties, one on chromosome 5 and one on
chromosome 6, are retained exactly. The tied markers have distinct
identifiers. No positional offset, marker renaming or deduplication is
applied to the physical map.

The three objects are stored together in \file{v-rapeseed.rda}
using xz compression at level 9. The STvcf \code{source},
\code{original.coding}, and \code{marker.selection} attributes record
provenance and the lossless conversion rule. These data are distinct
from the multi-environment height data loaded by \code{v-lmm-rapeseed}.
}
\value{
Loading the data creates \code{v.rapeseed.vcf},
\code{v.rapeseed.phe}, and \code{v.rapeseed.map}, as described
in the format section.
}
\source{
Data provided by Rod Snowdon (University of Giessen).
}
\keyword{datasets}
