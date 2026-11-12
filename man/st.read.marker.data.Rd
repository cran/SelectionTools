% copyright (C) 2001 Frisch/Maurer
%
\name{st.read.marker.data}
\alias{st.read.marker.data}
\title{st.read.marker.data()}
\description{
Reads marker data into a named SelectionTools marker data set and derives the marker statistics used by subsequent analyses.
}
\usage{
st.read.marker.data(filename, format, auxfiles = FALSE, data.set = "default")
}
\arguments{
\item{filename}{Character string naming the marker input file. Absolute paths are used directly; relative names are resolved using \code{st.input.dir}.}
   \item{format}{Character string selecting the input layout: \code{"l"} for list format, \code{"m"} for marker-major matrix format, \code{"t"} for transposed semicolon-separated matrix format, or \code{"n"} for NTSys/Plabsim allele-incidence format.}
   \item{auxfiles}{Logical. If \code{TRUE}, marker, individual, and genotype summary tables are returned directly from memory.}
   \item{data.set}{Character string naming the SelectionTools marker data set that is replaced after a successful read.}
}
%

\details{
Four input layouts are supported. Format \code{"l"} is a list or long format in which a record identifies an individual, a marker, and one allele. A homozygous genotype can therefore be represented by one allele record and a heterozygous genotype by two records for the same individual and marker. In this format a multi-digit token such as \code{12} is a single allele number 12.

Format \code{"m"} is a marker-major matrix. Individual names are read from the header and each subsequent row starts with a marker name followed by one genotype for each individual. Format \code{"t"} is the transposed, individual-major matrix and uses semicolon-separated fields. Format \code{"n"} is the NTSys/Plabsim allele-incidence representation in which rows have names of the form marker.allele followed by zero-one incidence values.

In matrix genotype fields, a single allele is interpreted as homozygous. A slash separates explicitly coded alleles, such as \code{12/14}. Compact two-character genotypes are interpreted as two one-character alleles, so \code{12} denotes genotype 1/2. Compact nucleotide genotypes are also recognized; the internal nucleotide codes are A=1, C=2, G=3, T=4, and D=5. Numeric biological allele codes must be integers from 1 through 254; they do not need to be consecutive within a marker. The tokens \code{-1}, \code{--}, \code{-}, \code{.}, and \code{NA} denote missing data. The \code{-1} spelling is retained for compatibility with historical SelectionTools marker files. If either allele of a diploid genotype is missing, the complete genotype is stored as missing.

The read is transactional. The new marker data, names, genotype matrix, and marker statistics are built and validated before replacing an existing data set. If the file cannot be read or is structurally invalid, the previous marker data set is left unchanged. Duplicate marker or individual identities that would make a matrix ambiguous are rejected. NTSys marker blocks must be contiguous and incidence values must be zero or one.

After a successful read the marker statistics are recomputed from the imported matrix. When \code{auxfiles=TRUE}, the individual, marker, and genotype summaries are returned directly from memory.
}

\value{
Invisibly returns a \code{list} when the requested auxiliary tables are produced, with components \code{individual.list}, \code{marker.list}, and \code{genotypes}, each a \code{data.frame}. Otherwise the return value is \code{NULL}; the main effect is loading marker data into \code{data.set}.
}

%
