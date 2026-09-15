% copyright (C) 2001 Frisch/Maurer
%
\name{st.marker.data.statistics}
\alias{st.marker.data.statistics}
\title{st.marker.data.statistics()}
\description{
Recomputes and optionally returns summary statistics for the current marker matrix.
}
\usage{
st.marker.data.statistics(data.set ="default", filename ="marker.stats",
                          ind = TRUE, mar = TRUE, gen = TRUE, auxfiles = TRUE)
}
\arguments{
\item{data.set}{Character string naming the marker data set.}
   \item{filename}{Character string retained for compatibility and used as part of the names of internal statistics files in the R session temporary directory. These temporary files are removed before the function returns.}
   \item{ind}{Logical. Request the individual missingness table.}
   \item{mar}{Logical. Request the marker statistics table.}
   \item{gen}{Logical. Request the genotype matrix table.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
}

\details{
Statistics are calculated directly from the two stored allele matrices. For each marker the function determines the distinct biological allele codes, allele counts, the proportion of missing allele copies, and expected heterozygosity calculated as one minus the sum of squared non-missing allele frequencies. The internal missing-allele sentinel is excluded from the allele-frequency denominator.

For each individual, the proportion of markers with a missing first allele is reported as the individual missingness statistic. These statistics are also the values used by marker-data restriction functions.

The marker summary distinguishes stored missingness from biological allele codes. Missing data are not counted as a biological allele for restrictions such as \code{NoAll.MAX}. Biological allele codes are in the range 1 through 254.

When requested, three files are produced: an individual summary, a marker summary, and a genotype matrix. The R wrapper reads these files back and removes the temporary files.
}
\value{
Invisibly returns a named \code{list} containing the requested marker-data summaries. Depending on \code{ind}, \code{mar}, and \code{gen}, components are \code{individual.list}, \code{marker.list}, and/or \code{genotypes}, each a \code{data.frame}; \code{NULL} is returned if no result can be produced.
}

%
%
