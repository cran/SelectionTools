\name{v-flintmaize2-vcf}
\alias{v-flintmaize2-vcf}
\alias{v.flintmaize2.vcf}
\docType{data}
\title{Mapped Common-Marker Data for Flint Maize Parents}
\description{STvcf marker data for 86 flint maize parents at 8072 mapped
markers that are polymorphic in both parental pools of the Melchinger and
Frisch (2023) hybrid-prediction example.}
\usage{data("v-flintmaize2-vcf", package = "SelectionTools")}
\format{A nine-component list of class \code{"STvcf"}:
\describe{
\item{\code{CHROM}}{8072 integer chromosome codes from 1 to 10.}
\item{\code{POS}}{8072 double-precision estimated genetic positions in cM.}
\item{\code{ID}}{8072 unique marker identifiers, retained from
\code{v.flintmaize.vcf}.}
\item{\code{REF}, \code{ALT}}{Allele labels \code{"1"} and \code{"2"}.}
\item{\code{individual}}{86 parental identifiers, unchanged from
\code{v.flintmaize.vcf} and matching \code{v.hybmaize2.phe}.}
\item{\code{allele1}, \code{allele2}}{8072 by 86 raw matrices. Raw 0
represents REF and raw 1 represents ALT. All calls are complete and
homozygous, so the matrices are identical.}
\item{\code{CHROM.levels}}{The chromosome labels \code{"1"}, \ldots,
\code{"10"}.}
}}
\details{The parental genotypes, allele coding and individual order are
retained from \code{v.flintmaize.vcf}. Markers are selected from the original
35478-marker panel by requiring both a map position and polymorphism among
the 86 flint parents and among the 123 dent parents. Both mapped parental
objects contain the same 8072 markers in the same chromosome/position order.
No additional marker filtering is needed to remove monomorphic loci.

The map matches \code{name_50k} to the existing marker identifiers with
\code{make.names}, uses \code{chr} as the chromosome, and retains
\code{pos_cM_estimated} as double-precision cM positions. The 10726 mapped
markers in the original panel contain 8072 that are polymorphic in both
pools. The supplied archive contributes the map; parental genotypes and
identifiers come from the existing SelectionTools data.

Two exact within-chromosome position ties are retained in the stored
objects. \code{st.load.vcf.data} separates ties internally by small
deterministic offsets and leaves the stored STvcf objects unchanged.
The \code{map.source} and \code{marker.selection} attributes describe the
map provenance and the marker-selection rule. The genotype source and
model reference are retained separately in the original attributes.}
\source{Parental genotypes: Technow et al. (2014), Genetics 197:1343--1355.
\doi{10.1534/genetics.114.165860}.
Model reference: Melchinger AE, Frisch M (2023). Genomic prediction in hybrid
breeding: II. Reciprocal recurrent genomic selection with full-sib and
half-sib families. Theoretical and Applied Genetics 136:203.
\doi{10.1007/s00122-023-04446-3}.
Map: supplied \file{maize2.tar.gz}, \file{map_HOH_pruned.csv}, columns
\code{name_50k}, \code{chr}, and \code{pos_cM_estimated}; these are exported
unchanged in \file{maize2.mmp}.}
\examples{
data("v-flintmaize2-vcf", package = "SelectionTools")
st.load.vcf.data(v.flintmaize2.vcf, data.set = "flint")
head(st.get.map(data.set = "flint"))
A.flint <- st.get.A("flint")
dim(A.flint)
}
\seealso{\code{\link{v-dentmaize2-vcf}}, \code{\link{v-hybmaize2-phe}},
\code{\link{v-flintmaize-vcf}}, \code{\link{st.load.vcf.data}},
\code{\link{st.set.simpop}}}
\keyword{datasets}
