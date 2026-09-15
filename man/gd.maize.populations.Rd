\name{gd.maize.populations}

\alias{gd.maize.populations}
\docType{data}

\title{SSR data from seven maize populations}

\description{
This data set consists of the SSR marker data of seven tropical maize populations.
From each populations 48 individuals were samples  analyzed with 85 polymorphic Single Sequence Repeat markers.
}

\usage{
data(gd.maize.populations)
}

\format{
The first five rows and columns look like:
   \tabular{cccccc}{
          \tab POL24.01 \tab POL24.02 \tab POL24.03 \tab POL24.04 \tab POL25.01\cr
nc130.139 \tab       0  \tab       0  \tab       1  \tab       1  \tab      NA\cr 
nc130.142 \tab       1  \tab       1  \tab       0  \tab       1  \tab       1\cr
nc130.145 \tab       0  \tab       1  \tab       0  \tab       0  \tab       0\cr
nc133.110 \tab       1  \tab       0  \tab       0  \tab       1  \tab       1\cr
nc133.148 \tab       0  \tab       0  \tab       0  \tab       0  \tab       0\cr
     }
Each column contains the marker genotype of an individual. The
identifier for an individual (the column name) contains the name of the
population to which the individual belongs and the name of the
individual, separated by a dot. For each allel of a marker there is a
row in the dataset, the rowname contains the name of the marker and the
name of the allele seperated by a dot.

The data matrix consists of zeros and ones, denoting whether a certain
allale was observed at an individual or not. Missing data are denoted
with NA.
}

\seealso{
  \code{\link{gd.genetic.distance}},
  \code{\link{gd.similarity.coefficient}},
  \code{\link{gd.maize.lines}},
  \code{\link{gd.maize.aflp}},
  \code{\link{gd.salad.aflp}}
}

\keyword{datasets}
