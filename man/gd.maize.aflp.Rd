\name{gd.maize.aflp}
\alias{gd.maize.aflp}

\alias{maize.aflp}
\docType{data}

\title{51 maize lines analyzed with 462 AFLP markers}

\description{
  This data set consists of AFLP data from 51 maize lines. 462 markers
  are used. This is dominant data. Each band is a marker, which has only
  one allele. Analysis of these dominant markers uses
  \code{gd.allow.zero.frequencies(1)}.
}

\usage{
data(gd.maize.aflp)
}

\format{
The first five rows and columns look like:
   \tabular{cccccc}{
              \tab l1.1 \tab  l2.1  \tab l3.1  \tab l4.1  \tab l5.1 \cr
E35M50+M001.1 \tab    1 \tab     1  \tab    1  \tab    0  \tab    1 \cr
E35M50+M002.1 \tab    0 \tab     0  \tab    0  \tab    0  \tab    1 \cr
E35M50+M003.1 \tab    0 \tab     1  \tab    0  \tab    0  \tab    0 \cr
E35M50+M004.1 \tab    1 \tab     1  \tab    1  \tab    0  \tab    1 \cr
E35M50+M005.1 \tab    0 \tab     0  \tab    1  \tab    0  \tab    1\cr
     }

Each column contains the marker genotype of an individual. The
identifier for an individual (the column name) contains the name of the
population to which the individual belongs and the name of the
individual, separated by a dot. For band there is a
row in the dataset, the rowname contains the name of the marker and 
 seperated by a dot a 1, indicating that each marker has only one allele. 


The data matrix consists of zeros and ones, denoting whether a certain
band was observed at an individual or not. Missing data are denoted
with NA.


For AFLPs the primer combination and the band is separated by a \code{+}.
}

\seealso{
  \code{\link{gd.genetic.distance}},
  \code{\link{gd.similarity.coefficient}},
  \code{\link{gd.maize.populations}},
  \code{\link{gd.maize.lines}},
  \code{\link{gd.salad.aflp}}
}

\keyword{datasets}
