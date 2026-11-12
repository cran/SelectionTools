\name{v-tropmaize-vcf}
\alias{v-tropmaize-vcf}
\alias{v.tropmaize.vcf}
\docType{data}
\title{Compacts STvcf Genotype and Genetic-Map Data for Tropical Maize}
\description{
A compact STvcf representation of genotypes and genetic linkage positions for
264 tropical maize inbred lines. It contains markers occurring in both the
tropical-maize marker matrix and the available genetic map.
}
\usage{
data("v-tropmaize-vcf")
}
\format{
A list of class \code{"STvcf"} with 1076 marker rows and 264 individuals:
\describe{
\item{\code{CHROM}}{Integer running chromosome codes 1 through 10.}
\item{\code{POS}}{Double-precision genetic positions in cM retained from the
available linkage map without intentional rounding.}
\item{\code{ID}}{Character vector containing 1076 unique marker IDs.}
\item{\code{REF}}{Character vector containing \code{"1"} for every marker.}
\item{\code{ALT}}{Character vector containing \code{"2"} for every marker.}
\item{\code{individual}}{Character vector \code{"1"}, ..., \code{"264"}.}
\item{\code{allele1}}{A 1076 by 264 raw matrix. Raw 0 denotes REF, raw 1
denotes ALT, and raw 255 denotes missing.}
\item{\code{allele2}}{A 1076 by 264 raw matrix with the same coding.}
\item{\code{CHROM.levels}}{Character vector \code{"1"}, ..., \code{"10"}.}
}
}
\details{
The original tropical-maize genotype object contains 1135 genotyped markers
and the available linkage map contains 7140 entries. This STvcf object
contains their 1076-marker intersection. Genotype-only markers lacking a map
position and map-only loci lacking genotype data are not included.

POS contains the original double-precision cM positions from the distributed
linkage map. The positions are not rounded to integer cM. Exact ties are
permitted in STvcf; when installed with \code{st.load.vcf.data()}, ties are
separated internally by small deterministic offsets because the existing
SelectionTools map and simulation structures require distinct ordered locus
positions.

REF label \code{"1"} is VCF allele index 0 and ALT label \code{"2"} is VCF
allele index 1. Loading maps these indices to internal SelectionTools alleles
1 and 2 and maps raw 255 to the internal SelectionTools missing representation.

The marker and individual identifiers obey the STvcf/SelectionTools name
convention used by the converters: identifiers are non-empty and contain no
whitespace or control characters. Ordinary punctuation is permitted. The
tropical-maize marker IDs happen to be alphanumeric and the individual IDs are
the digit strings 1 through 264.

Marker statistics and individual missing-data frequencies after loading are
calculated from the 1076 markers in this object. Applying NoAll.MAX=2,
MaMis.MAX=0.1, ExHet.MIN=0.1, and InMis.MAX=0.1 leaves 828 markers and 230
individuals, the intended tropical-maize analysis data set.
}
\source{
Jose Crossa, Gustavo de los Campos, Paulino Perez, Daniel Gianola, 
Juan Burgueno, Jose Luis Araus, Dan Makumbi, Ravi P Singh, Susanne Dreisigacker, 
Jianbing Yan, Vivi Arief, Marianne Banziger, Hans-Joachim Braun, 
Prediction of Genetic Values of Quantitative Traits in Plant Breeding Using Pedigree and Molecular Markers, 
Genetics, Volume 186, Issue 2, 1 October 2010, Pages 713-724, 
\doi{10.1534/genetics.110.118521}
}
\keyword{datasets}
