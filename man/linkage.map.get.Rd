\name{linkage.map.get}
\alias{linkage.map.get}
\alias{population.marker.map}
\alias{extend.marker.map}
\title{Returns a linkage map used by the simulation routines}

\description{
Returns the currently active simulation linkage map.
}

\usage{
linkage.map.get()

population.marker.map(n.markers)

extend.marker.map(n.markers)
}

\details{
The result contains the map points used by the simulation backend, including chromosome, position, locus name, and class, in the current simulation map order.
}

\value{
A \code{data.frame} with columns \code{chrom}, \code{pos}, \code{name}, and \code{class}, describing the active simulation linkage map.
}

\arguments{\item{n.markers}{Total number of map points after the operation. For replacement, at least the number of chromosomes in the defined genome; for extension, at least the current map size.}}

\section{Random simulation marker map}{
Replaces the simulation marker map with a specified number of random marker positions.

Places at least one marker on each chromosome and allocates remaining markers with probabilities proportional to chromosome length. Positions are uniform within chromosomes and sorted physically. Marker names are \code{M1}, \code{M2}, and so on. Uses the R random-number generator; \code{set.seed} provides reproducibility. The chromosome segments and stored pedigrees of existing simulation populations are retained. Genotypes at the new positions are obtained when populations are exported to marker datasets.

The active linkage map, invisibly.

}

\section{Extension of an existing marker map}{
\code{extend.marker.map} retains every existing locus with its identifier,
position and class, and adds random marker positions until the map contains
\code{n.markers} points in total. New positions are uniform over genome length;
new identifiers avoid existing names. The result is sorted physically within
chromosomes. Chromosome lengths, inherited segments and stored pedigrees are
preserved. Genotypes at the new loci are derived from those inherited segments,
not additional observed founder genotypes. The operation invalidates cached
simulation marker genotypes, which are rebuilt when requested. As with map
replacement, simulation effect-map definitions and cached evaluated values are
cleared; this operation does not change marker-dataset model results.

An equal count is a no-op and does not advance the random-number generator.
A smaller count is rejected without changing the map. Uses the R random-number
generator; its seed determines reproducibility. Returns the active map invisibly.
}
