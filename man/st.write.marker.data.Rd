% copyright (C) 2001 Frisch/Maurer
%
\name{st.write.marker.data}
\alias{st.write.marker.data}
\title{st.write.marker.data()}

\description{
Writes marker data, and where applicable map and genome information, in one of the supported external representations.
}
\usage{
st.write.marker.data(format = "m", lfilename = NULL, mfilename = NULL,
                     nfilename = NULL, ifilename = "", f.ind = 0,
                     l.ind = 0, add.inum = ".1", auxfiles = FALSE,
                     data.set = "default")
}
\arguments{
\item{format}{Character string selecting \code{"l"}, \code{"m"}, or \code{"n"}. The C writer also supports \code{"a"} as append-list mode.}
   \item{lfilename}{Base name for list-format output. Required when \code{format = "l"} or \code{format = "a"}.}
   \item{mfilename}{Base name for matrix-format marker and map output. Required when \code{format = "m"}.}
   \item{nfilename}{Base name for NTSys/Plabsim marker, map, and genome-parameter output. Required when \code{format = "n"}.}
   \item{ifilename}{Optional file containing individual names to include in the output.}
   \item{f.ind}{First individual position to write when selecting by range. Positions are interpreted in the current data-set order.}
   \item{l.ind}{Last individual position to write when selecting by range.}
   \item{add.inum}{Character string appended to individual names in NTSys/Plabsim output.}
   \item{auxfiles}{Logical flag retained by the interface. The marker writer itself writes the files implied by \code{format}.}
   \item{data.set}{Character string naming the marker data set to write.}
}

\details{
Format \code{"l"} writes list records consisting of individual, marker, and allele. A heterozygous genotype is written as two allele records. The internal append variant \code{"a"} appends list records to the list output file without writing a new header.

Format \code{"m"} writes a marker-major genotype matrix. When a linkage map is available, a corresponding map file is also written. Genotypes are written with an explicit slash between the two allele codes, which preserves multi-digit allele numbers unambiguously.

Format \code{"n"} writes the Plabsim/NTSys group of files: allele-incidence marker data, map data in Morgan units, and chromosome-length information. Individual names are adapted to the restrictions of this representation and \code{add.inum} is appended to the written name. Marker names containing a period are not suitable for this output because the period separates marker name and allele number. Name transformations that would make two output individuals indistinguishable are rejected.

The writer does not modify the in-memory marker data set. In particular, writing a data set without a linkage map uses the marker matrix directly and does not create a synthetic map or alter map dimensions. Unsupported format codes are rejected.

Individuals can be restricted either by a positional range or by a file of individual names. These output selections affect only the written files and do not alter the stored marker data.
}
\value{
\code{NULL}, returned invisibly. The function is called for its side effect of writing marker-data files in the selected output format.
}

%
%
