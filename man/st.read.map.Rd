% copyright (C) 2001 Frisch/Maurer
%
\name{st.read.map}
\alias{st.read.map}
\title{st.read.map()}
\description{
Reads a linkage map for an existing marker data set, matches loci by marker name, and restricts the marker matrix to markers represented in both sources.
}
\usage{
st.read.map(filename, skip = 0, format = "cpms", m.stretch = 1,
            auxfiles = TRUE, data.set = "default")
}
\arguments{
\item{filename}{Character string naming the map input file. Absolute paths are used directly; relative names are resolved using \code{st.input.dir}.}
   \item{skip}{Non-negative integer giving the number of initial input lines to ignore.}
   \item{format}{Character string specifying \code{"cpms"} or \code{"mcp"}.}
   \item{m.stretch}{Finite numeric multiplier applied to every map position after reading. A value of 1 retains positions already in the desired map unit.}
   \item{auxfiles}{Logical. If \code{TRUE}, the matched and sorted map is returned directly from memory.}
   \item{data.set}{Character string naming the marker data set to which the map is attached.}
}
%

\details{
A marker data set must already be loaded. The supported map layouts are \code{"cpms"}, containing chromosome, position, marker name, and class, and \code{"mcp"}, containing marker name, chromosome, and position. In \code{"mcp"} input the class is assigned internally. The first \code{skip} lines can be ignored before map parsing. Map positions are multiplied by \code{m.stretch}.

Marker matching is by marker name. Loci present in the map but absent from the marker matrix are silently ignored. Marker rows present in the marker matrix but absent from the map are silently removed. This behavior is intentional: after a successful read the marker data set represents the intersection of marker matrix and linkage map. The remaining genotype rows are reordered into linkage-map order and marker statistics are recalculated. Phenotypic data associated with individuals are preserved.

Duplicate marker names in either the linkage map or marker matrix are rejected because they make the mapping ambiguous. Equal positions on the same chromosome are separated deterministically while preserving their relative order and without moving a tied locus beyond the next genuinely distinct map position.

The operation is transactional. Parsing, matching, construction of the selected genotype matrix, linkage-map structures, and marker statistics are completed before the existing marker data are replaced. On failure, the previous marker matrix, map, phenotype vector, and derived state remain available. On successful replacement, model- and marker-order-dependent derived structures belonging to the old marker set are invalidated.
}

\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{Pos}, \code{Name}, and \code{Class} when \code{auxfiles = TRUE} and map import succeeds; otherwise \code{NULL}. The main effect is updating the map associated with \code{data.set}.
}

%
