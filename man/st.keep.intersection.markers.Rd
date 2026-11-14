% copyright (C) 2001 Frisch/Maurer
%
\name{st.keep.intersection.markers}
\alias{st.keep.intersection.markers}
\title{Keeps the Marker Intersection of Two SelectionTools Data Sets}
\description{
Restricts two SelectionTools marker data sets to markers present in both data sets.
}
\usage{
st.keep.intersection.markers(data.set1, data.set2)
}
\arguments{
\item{data.set1}{One non-empty character string naming the first marker data set.}
\item{data.set2}{One non-empty character string naming the second marker data set.}
}
\details{
Marker matching is by exact marker name. Both data sets retain only the intersection, and the second data set is reordered so that the common markers occur in the same order as in \code{data.set1}. Duplicate marker names are rejected.

Individual names and genotype calls are retained. Phenotypic data are retained when present. If a data set has a linkage map, the corresponding map entries for the retained markers are preserved and reordered with the marker matrix. Marker statistics are recalculated. Marker-order-dependent derived structures such as fitted model objects, haplotype definitions, and cached chromosome statistics are invalidated when a data set must be rebuilt.

The operation is transactional across the two data sets: both subset replacements are constructed before either existing data set is changed. If construction fails, both original data sets remain unchanged. If the two arguments name the same data set, no change is made.
}
\value{
Invisibly returns the number of markers in the intersection.
}
