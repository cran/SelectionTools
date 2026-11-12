\name{st.data.combine}
\alias{st.data.combine}
\title{Combine Individuals from Marker Data Sets}
\description{
Creates a marker data set by appending individuals from source data sets
having exactly the same marker identifiers. Validation, matching, copying
and storage are implemented in C.
}
\usage{
st.data.combine(target.data.set, source.data.sets)
}
\arguments{
\item{target.data.set}{One nonempty character string naming the target
marker data set. It must differ from every source name.}
\item{source.data.sets}{A nonempty list of single dataset-name strings,
or a character vector of dataset names. All sources must already contain
valid marker data.}
}
\details{
Every source must have exactly the same set of unique marker identifiers.
Marker order may differ: matching uses marker names, and the output marker
order follows the first source. No intersection, union, automatic dropping
of markers or allele recoding is performed. The observed allele sets may
differ between sources, including fixed alleles and missing observations.
Allele codes must already be compatible. For nucleotide VCF inputs,
\code{strictACTG=TRUE} in \code{read.vcf} or \code{st.load.vcf.data} provides
fixed codes independent of REF/ALT order.

Individuals appear in source-list order and retain their order within each
source. Individual identifiers must be unique across the combined data;
duplicate identifiers, repeated source names and missing source names cause
an error. Names and pedigrees are copied independently of the sources.
Genotype order and missing values are preserved. The output is marked as
phased only if every source is marked as phased.

The linkage map is copied from the first source, including its map order
and mapping to marker rows. Maps from other sources are not combined or
compared. If the first source has no linkage map, neither does the output.

Available phenotypes are appended by individual. If only some sources
contain phenotypes, individuals from the other sources receive missing
phenotypes. Marker statistics are recalculated, and the usual automatic
packed-storage policy applies. Fitted effects, prediction state, relationship
matrices and derived workspaces are not inherited.

Source data sets are unchanged. A pre-existing target is replaced only after
the combined data have been prepared successfully. A target referenced by
a hybrid data set cannot be replaced. Invalid input leaves an existing
target unchanged and does not create a new empty target.

Large genotype copies use the thread budget set by \code{st.set.num.threads}
when the package has OpenMP support.
}
\value{
Integer zero, invisibly, on success. Invalid or incompatible inputs cause
an error. The combined marker data set is stored as a side effect.
}
