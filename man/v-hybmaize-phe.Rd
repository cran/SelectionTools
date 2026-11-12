\name{v-hybmaize-phe}
\alias{v-hybmaize-phe}
\alias{v.hybmaize.phe}
\docType{data}
\title{Phenotype Data for Maize Hybrids}
\description{Phenotype records for 1254 crosses between dent and flint maize parents
in the Melchinger and Frisch (2023) hybrid-prediction example data set. Parental identifiers match the accompanying marker data.}
\usage{data("v-hybmaize-phe")}
\format{A data frame with 1254 rows and six columns:
\describe{
\item{\code{hybrid}}{Hybrid identifier, with parents separated by a period.}
\item{\code{dent}}{Dent parent identifier.}
\item{\code{flint}}{Flint parent identifier.}
\item{\code{GY}}{Grain yield phenotype, as supplied in the source data.}
\item{\code{GM}}{Grain moisture phenotype, as supplied in the source data.}
\item{\code{hy}}{Hybrid identifier, with parents separated by a colon.}
}}
\details{The \code{dataset.name}, \code{model.reference} and \code{source} attributes retain the example label, model reference and original data provenance separately.}
\source{Model reference: Melchinger AE, Frisch M (2023). Genomic prediction in hybrid breeding: II. Reciprocal recurrent genomic selection with full-sib and half-sib families. Theoretical and Applied Genetics 136:203. \doi{10.1007/s00122-023-04446-3}.
Original parental markers and hybrid records: Technow et al. (2014), Genetics 197:1343--1355. \doi{10.1534/genetics.114.165860}.}
\keyword{datasets}
