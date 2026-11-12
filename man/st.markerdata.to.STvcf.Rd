\name{st.markerdata.to.STvcf}
\alias{st.markerdata.to.STvcf}
\title{Converts a SelectionTools Marker Data Set to STvcf}
\description{
Converts marker genotypes and the linkage map of an existing SelectionTools
marker data set to the compact diploid STvcf representation.
}
\usage{
st.markerdata.to.STvcf(data.set="default")
}
\arguments{
\item{data.set}{Character string naming the SelectionTools marker data set to
be converted.}
}
\details{
A stored physical map is retained in the exported \code{physical.map}
attribute. A dataset with a physical map but no genetic map is exported
with missing \code{CHROM}/\code{POS} and empty \code{CHROM.levels};
base-pair positions are never substituted for genetic positions in cM.

The optional \code{allele.codes} attribute is a list of integer coding vectors, named and ordered by marker ID. It records the internal allele codes corresponding to REF and successive ALT entries. STvcf exchange retains this attribute so fixed or subsetted loci keep the same coding as their fitted marker reference. The loader validates it before replacing a destination. Conversions preserve it; a manually edited marker list must keep it aligned. It is not part of plain-text VCF exchange.

This function supplies the marker-data to STvcf step of the reverse data
chain. With the existing \code{st.get.simpop()} function, a simulation
population can first be returned to the ordinary SelectionTools marker-data
structure and can then be converted to STvcf without introducing a separate
simulation-to-STvcf bridge.

The marker data set must contain genotype data and either a usable
genetic linkage map or a stored physical map. When a genetic map is
present, every marker-matrix column must have exactly one active genetic
map location. Map-only
locations are ignored. Missing or duplicate marker locations, inconsistent
marker and map names, invalid chromosome numbers, non-finite positions, and
negative positions cause conversion to fail.

The current marker order is retained. Map positions are returned as the
current double-precision cM values without intentional rounding. Current
numeric chromosome numbers are converted to consecutive STvcf running codes
1, 2, ..., in order of first appearance, while \code{CHROM.levels} records
the corresponding numeric chromosome labels as character strings.

Current marker and individual names become STvcf identifiers. They must be
unique and non-empty. Marker names must contain fewer than 256 bytes; individual names use owned dynamic strings subject to the R string-size limit. Whitespace and control
characters are rejected; ordinary punctuation such as dots, underscores, and
hyphens is permitted. Names are never silently modified during this conversion.

Simulation exports retain each individual's stable name and complete stored pedigree. Imported founders keep their imported names. New breeding individuals use their creation population and serial number followed by the stored pedigree expression. Copying or renaming a population does not regenerate these identifiers.

The SelectionTools marker structure stores integer allele states but not the
original textual VCF REF/ALT labels. The converter therefore normalizes the
observed non-missing allele states separately for each marker. States are
ordered by their current integer values. The first state becomes
\code{REF="1"}; further states become \code{ALT="2,3,..."}; genotypes are
converted to VCF indices 0, 1, 2, ... in that order.

For marker data originally loaded from STvcf, internal alleles are already 1
for REF, 2 for ALT1, 3 for ALT2, and so on. Transfer to and from simulation
preserves those integer states. For marker data created by other routes,
distinct genotype states are retained but historical integer labels are
normalized.

An all-missing marker is retained with \code{REF="1"}, \code{ALT="."}, and
missing genotypes. A marker with one observed non-missing state is likewise
written with \code{ALT="."} and uses VCF allele index 0. If either homologue
of a marker-data genotype is missing, the exported STvcf genotype is
completely missing, consistent with the STvcf partial-missing convention.

The two current marker-data allele arrays are retained in order as
\code{allele1} and \code{allele2}. The returned STvcf object receives a
logical \code{phased} attribute from the compiled marker data set. For data
obtained through \code{st.get.simpop()}, the two arrays correspond to the two
simulation homologues and the attribute is \code{TRUE}. No new phasing
calculation is performed.

The current SelectionTools marker representation permits at most 254
biological allele states per marker. After normalization their VCF indices are
therefore at most 253, so this converter writes raw STvcf allele matrices with
raw value 255 for missing. Marker genotypes that violate the SelectionTools
allele-range or missingness invariants are rejected.

The linkage-map class is not part of STvcf. Exact map-position ties may have
already been dispersed by the SelectionTools linkage-map or simulation
machinery; this function returns the current stored positions and cannot
reconstruct an earlier duplicate position that is no longer represented.

The optional \code{pedigree} attribute contains the versioned \code{name}/\code{pedigree} data frame. When an active symbolic-A fit is stored, \code{prediction.state} retains its versioned fitted snapshot and source role for STvcf/RDS persistence. These attributes leave the nine-component STvcf layout unchanged. Standalone observed responses, legacy effects and selection values, and simulation GValue/PValue arrays are not exported as dataset fields. The fitted snapshot contains the response and effect information needed for prediction; refitting after import requires loading the desired response into the marker dataset.
}
\value{
A \code{list} of class \code{"STvcf"} containing the current marker data, individual identifiers, marker metadata and chromosome-label mapping. If conversion fails, \code{NULL} is returned.
}
