\name{v-flintmaize-vcf}
\alias{v-flintmaize-vcf}
\alias{v.flintmaize.vcf}
\docType{data}
\title{STvcf Marker Data for Flint Maize Parents}
\description{Lossless STvcf marker data for 86 parents and 35478 markers in
the Melchinger and Frisch (2023) hybrid-prediction example data set.}
\usage{data("v-flintmaize-vcf", package = "SelectionTools")}
\format{A nine-component list of class \code{"STvcf"}:
\describe{
\item{\code{CHROM}}{35478 missing integer chromosome codes.}
\item{\code{POS}}{35478 missing double-precision map positions.}
\item{\code{ID}}{Unique marker identifiers in original marker order.}
\item{\code{REF}, \code{ALT}}{Allele labels \code{"1"} and \code{"2"}.}
\item{\code{individual}}{86 parental identifiers, matching the hybrid records.}
\item{\code{allele1}, \code{allele2}}{35478 by 86 raw matrices. Raw 0 denotes the reference allele and raw 1 the alternative allele. The two matrices agree because the supplied parents are homozygous.}
\item{\code{CHROM.levels}}{Empty character vector, explicitly indicating marker-only data without a linkage map.}
}}
\details{\code{st.load.vcf.data} loads \code{v.flintmaize.vcf} as a SelectionTools marker data set. The supplied object preserves the original marker genotypes and identifiers. No linkage map is supplied and no artificial positions are assigned. Relationship matrices depend on the marker genotypes and do not require a map; map-dependent analyses require separately supplied genuine map information.

For \code{st.get.A}, \code{st.restrict.marker.data(ExHet.MIN = 1e-8, NoAll.MAX = 2, data.set = "flintmaize")} removes markers that are monomorphic within the flint pool. This retains 28240 polymorphic biallelic markers. The relationship constructor reuses the existing marker-to-Z and VanRaden A-matrix functions. The R object \code{v.flintmaize.vcf} retains all 35478 markers.

The \code{dataset.name}, \code{model.reference} and \code{source} attributes distinguish the Melchinger and Frisch (2023) example label from the original marker-data provenance.}
\source{Model reference: Melchinger AE, Frisch M (2023). Genomic prediction in hybrid breeding: II. Reciprocal recurrent genomic selection with full-sib and half-sib families. Theoretical and Applied Genetics 136:203. \doi{10.1007/s00122-023-04446-3}.
Original parental markers and hybrid records: Technow et al. (2014), Genetics 197:1343--1355. \doi{10.1534/genetics.114.165860}.}
\examples{
data("v-flintmaize-vcf", package = "SelectionTools")
st.load.vcf.data(v.flintmaize.vcf, data.set = "flintmaize")
st.restrict.marker.data(ExHet.MIN = 1e-8, NoAll.MAX = 2,
                        data.set = "flintmaize")
Aflint <- st.get.A("flintmaize")
dim(Aflint)
}
\seealso{\code{\link{v-dentmaize-vcf}}, \code{\link{v-hybmaize-phe}}, \code{\link{st.load.vcf.data}}, \code{\link{st.restrict.marker.data}}, \code{\link{st.get.A}}}
\keyword{datasets}
