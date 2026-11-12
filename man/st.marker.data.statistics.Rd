% copyright (C) 2001 Frisch/Maurer
%
\name{st.marker.data.statistics}
\alias{st.marker.data.statistics}
\title{Summarizes a SelectionTools Marker Data Set}
\description{
Recomputes and optionally returns summary statistics for the current marker matrix.
}
\usage{
st.marker.data.statistics(data.set ="default", filename ="marker.stats",
                          ind = TRUE, mar = TRUE, gen = TRUE, auxfiles = TRUE)
}
\arguments{
\item{data.set}{Character string naming the marker data set.}
   \item{filename}{Character string retained for compatibility. The statistics are returned directly from memory.}
   \item{ind}{Logical. TRUE requests the individual missingness table.}
   \item{mar}{Logical. TRUE requests the marker statistics table.}
   \item{gen}{Logical. TRUE requests the genotype matrix table.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
}

\details{
Statistics are calculated directly from the two stored allele matrices. 
For each marker the function determines the distinct biological allele codes, allele counts, 
the proportion of missing allele copies, and expected heterozygosity calculated as one minus 
the sum of squared non-missing allele frequencies.

For each individual, the proportion of markers with a missing first allele is reported as the individual missingness statistic. 
These statistics are also the values used by marker-data restriction functions.

The returned \code{individual.list} has columns \code{Name} and \code{InMis},
where \code{InMis} is the individual missing-marker frequency.

The returned \code{marker.list} contains \code{Name}, \code{NoAll},
\code{MaMis}, and \code{ExHet}. \code{NoAll} is the number of observed
non-missing biological alleles, \code{MaMis} is the marker missing-allele
frequency, and \code{ExHet} is the expected heterozygosity. The remaining
columns contain allele copy counts. If missing alleles occur anywhere in the
data set, column \code{AM} gives the number of missing allele copies at each
marker. Columns \code{A1}, \code{A2}, and so on give the counts of the
corresponding biological allele codes observed in the data set. A marker for
which a particular code is absent has count zero in that column.

The returned \code{genotypes} component is the complete marker-by-individual
genotype table. Its first column identifies the marker and the remaining
columns are named by individual. Diploid entries are printed as \code{allele1/allele2}. If the first
stored allele is missing, the entry is printed as \code{--}. 
\code{1} represents the REF allele, all ALT alleles are numbered starting with 2.
}
\value{
Invisibly returns a named \code{list} containing the requested marker-data summaries. Depending on \code{ind}, \code{mar}, and \code{gen}, components are \code{individual.list}, \code{marker.list}, and/or \code{genotypes}, each a \code{data.frame}; \code{NULL} is returned if no result can be produced.
}

\examples{
data("v-tropmaize-vcf", package = "SelectionTools")
marker.frame <- st.STvcf.to.dataframe(v.tropmaize.vcf)[1:12, c(1:6, 7:14)]
vcf <- st.dataframe.to.STvcf(marker.frame)
st.load.vcf.data(vcf, data.set = "maize.example")
statistics <- st.marker.data.statistics(data.set = "maize.example")
statistics$genotypes
}

\seealso{
\code{\link{st.restrict.marker.data}}, \code{\link{st.plot.gene.diversity}}
}

%
%
