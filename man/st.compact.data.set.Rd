% copyright (C) 2001 Frisch/Maurer
%
\name{st.compact.data.set}
\alias{st.compact.data.set}
\title{Compact a SelectionTools Marker Data Set}
\description{
Releases calculated results and reusable workspace allocations associated with
one SelectionTools marker data set while retaining its marker data and linkage
map.
}
\usage{
st.compact.data.set(data.set = "default")
}
\arguments{
  \item{data.set}{Character string naming the marker data set to compact.}
}
\details{
The retained core consists only of individual and marker names, the two allele
arrays containing the marker genotypes, the data-set phase state, and the
linkage map including its marker-to-data indices.

All other dynamically allocated data-set state is discarded. This includes
marker statistics, phenotypic and genotypic values, effects and model results,
design and relationship matrices, haplotype-block and cross structures,
chromosome-specific recombination/LD workspaces, and reusable temporary
buffers. Consequently columns such as \code{effects}, \code{phenotypic.value}, and \code{genotypic.value} in
\code{st.data.set.list()} can change from \code{TRUE} to \code{FALSE} after
compaction.

Values that are required again must be recalculated or reloaded by the
corresponding SelectionTools function.
}
\value{
Invisibly returns the estimated amount of work memory released, in MiB. The
estimate excludes memory allocator overhead.
}
