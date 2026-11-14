\name{v-hybmaize2-phe}
\alias{v-hybmaize2-phe}
\alias{v.hybmaize2.phe}
\docType{data}
\title{Hybrid Phenotypes for the Mapped Maize Parental Data}
\description{Phenotype records for 1254 crosses between the 86 flint and
123 dent parents in \code{v.flintmaize2.vcf} and \code{v.dentmaize2.vcf}.}
\usage{data("v-hybmaize2-phe", package = "SelectionTools")}
\format{A data frame with 1254 rows and six columns:
\describe{
\item{\code{hybrid}}{Hybrid identifier, with parents separated by a period.}
\item{\code{dent}}{Dent parent identifier matching \code{v.dentmaize2.vcf}.}
\item{\code{flint}}{Flint parent identifier matching \code{v.flintmaize2.vcf}.}
\item{\code{GY}}{Grain yield phenotype, as supplied in the source data.}
\item{\code{GM}}{Grain moisture phenotype, as supplied in the source data.}
\item{\code{hy}}{Hybrid identifier, with parents separated by a colon.}
}}
\details{This object is an unchanged copy of \code{v.hybmaize.phe}, including
all values, identifiers, row order, column types and attributes. Every parent
has marker data in the corresponding mapped parental object. The marker
subset changes the genomic information available for fitting, so a model
using the mapped panel is fitted to that panel rather than borrowing
estimates from the larger marker panel.

Direct storage of these records requires marker-data names
\code{"flint"} and \code{"dent"}, matching the parental column names in the
stored performance frame.}
\source{Melchinger AE, Frisch M (2023). Genomic prediction in hybrid
breeding: II. Reciprocal recurrent genomic selection with full-sib and
half-sib families. Theoretical and Applied Genetics 136:203.
\doi{10.1007/s00122-023-04446-3}.
Original parental markers and hybrid records: Technow et al. (2014),
Genetics 197:1343--1355. \doi{10.1534/genetics.114.165860}.}
\seealso{\code{\link{v-flintmaize2-vcf}}, \code{\link{v-dentmaize2-vcf}},
\code{\link{v-hybmaize-phe}}, \code{\link{st.hybrid.data.set}},
\code{\link{st.lmm}}}
\keyword{datasets}
