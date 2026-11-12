% copyright (C) 2001 Frisch/Maurer
%
\name{st.read.performance.data}
\alias{st.read.performance.data}
\title{st.read.performance.data()}
\description{
Reads phenotypic values for a marker data set and retains marker-data individuals with a matching phenotype record, including explicit missing values.
}
\usage{
st.read.performance.data(filename, out.filename = "y.vector", auxfiles = TRUE,
                         data.set ="default")
}
\arguments{
\item{filename}{Character string naming the phenotype input file. Relative names are resolved using \code{st.input.dir}.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Character string naming the marker data set whose individuals are matched to the phenotype file.}
}

\details{
The input contains individual identifiers and phenotypic values. The marker matrix must already be loaded because individual names are matched against the marker-data individuals.

Phenotype records for individuals that are not present in the marker matrix are silently ignored. This is a normal situation and does not generate a warning. Conversely, individuals present in the marker matrix for which no phenotypic record is found are removed from the marker data set. Thus a successful call leaves the individual intersection of marker and phenotype data. Marker order and the linkage map are preserved, while marker statistics are recalculated for the retained individuals.

If a known individual occurs more than once in the phenotype input, the last record is used. A phenotype for a known individual must be a valid finite numeric value or the token \code{NA}. An explicit \code{NA} record keeps the individual and its marker genotypes, with a missing phenotype; an individual without a matching record is removed. \code{NaN}, infinities and other nonnumeric tokens are rejected. Unknown individuals are ignored before their value field is interpreted. If no phenotype record matches any marker-data individual, the operation fails rather than creating an empty marker data set.

The operation is transactional. The reduced marker matrix, names, phenotype vector, map copy, and marker statistics are prepared before replacing the existing data set. A read or allocation failure therefore leaves the previous marker data set unchanged.
}

\value{
Invisibly returns a \code{data.frame} with columns \code{i} (individual identifier) and \code{y} (performance value) when \code{auxfiles = TRUE} and the operation succeeds; otherwise \code{NULL}. The main effect is updating performance data in \code{data.set}.
}

%
%
