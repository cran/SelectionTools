\name{st.read.marker.data.df}
\alias{st.read.marker.data.df}
\title{Load marker data from an R data frame or matrix}
\description{
Loads diploid marker genotypes into a named SelectionTools data set.
}
\usage{
st.read.marker.data.df(markerdata, format = NULL, data.set = "default")
}
\arguments{
\item{markerdata}{For \code{format="list"}, a data frame whose first
three columns are character vectors containing individual names, marker
names, and allele codes. For \code{format="matrix"}, a numeric or character
matrix or data frame with individuals in rows and markers in columns.
Explicit row names identify individuals and column names identify markers.
Both sets of names must be unique, nonempty, and contain no whitespace.}
\item{format}{Either \code{"list"} or \code{"matrix"}.
The default, \code{NULL}, selects list format.}
\item{data.set}{Character string naming the data set to create or replace.}
}
\details{
In list format, each row supplies an allele record for one individual and
marker. Alleles may be supplied as separate records or as pairs separated
by a slash, for example \code{"1/2"}. Biological allele codes range from 1 to 254. Missing values and the
tokens \code{"-1"}, \code{"--"}, \code{"-"}, \code{"."}, and \code{"NA"}
represent missing alleles. If either allele is missing, the complete
genotype is missing. Individual and marker names are sorted in list format.

In matrix format, rows represent individuals and columns represent markers.
Row names identify individuals and column names identify markers.

Genotypes may be supplied as dosages: \code{0} denotes a homozygous
major-allele genotype, \code{1} a heterozygote, and \code{2} a homozygous
minor-allele genotype. These correspond to allele pairs \code{1/1},
\code{1/2}, and \code{2/2}, respectively.

Alternatively, genotypes may be supplied as character allele pairs:
\code{"1/1"}, \code{"1/2"}, \code{"2/1"}, or \code{"2/2"}.
Use a pipe to indicate known phase, for example \code{"1|2"}.
The supplied allele order is preserved for both separators.

The data set is marked as phased only when all non-missing genotypes
use pipe-separated allele pairs. Mixing these with slash-separated
pairs or dosages marks the data set as unphased. Missing genotypes
do not affect this classification; entirely missing data are unphased.

Allele codes 1 and 2 refer to the user-designated major and minor
alleles; these labels are not inferred from observed frequencies.
R missing values, \code{"--"}, and \code{"-"} represent missing
genotypes. Partially missing pairs and other allele codes are not
supported. Matrix format preserves the supplied individual and marker order.
Invalid matrix input leaves an existing data set unchanged.
}
\value{
Returns \code{NULL} invisibly. The function loads marker data into
\code{data.set}.
}
