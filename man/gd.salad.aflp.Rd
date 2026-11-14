\name{gd.salad.aflp}

\alias{gd.salad.aflp}
\docType{data}

\title{Jasminas 44 salads analyzed with 108 AFLP markers}

\description{
  This data set consists of AFLP data from 44 salad lines. 108 markers
  are used. This is dominant data. Each band is a marker, which has only
  one allele. Analysis of these dominant markers uses
  \code{gd.allow.zero.frequencies(1)}.
}

\usage{
data(gd.salad.aflp)
}

\format{
The first five rows and columns look like:
   \tabular{cccccc}{
        \tab l34.1 \tab  l44.1 \tab l45.1 \tab l38.1 \tab l33.1\cr
p1+6.1  \tab    1  \tab     0  \tab    0  \tab    1  \tab    1 \cr
p2+1.1  \tab    1  \tab     1  \tab    1  \tab    1  \tab    1 \cr
p2+6.1  \tab    1  \tab     1  \tab    1  \tab    1  \tab    1 \cr
p2+7.1  \tab    1  \tab     1  \tab    1  \tab    1  \tab    1 \cr
p2+8.1  \tab    1  \tab     1  \tab    1  \tab    1  \tab    1 \cr
     }

Each column contains the marker genotype of an individual. The
identifier for an individual (the column name) contains the name of the
population to which the individual belongs and the name of the
individual, separated by a dot. For band there is a
row in the dataset, the rowname contains the name of the marker and 
 seperated by a dot a 1, indicating that each marker has only one allele. 


The data matrix consists of zeroes and ones, denoting whether a certain
band was observed at an individual or not. Missing data are denoted
with NA.


For AFLPs the primer combination and the band is separated by a \code{+}.

}

\seealso{
  \code{\link{gd.genetic.distance}},
  \code{\link{gd.similarity.coefficient}},
  \code{\link{gd.maize.populations}},
  \code{\link{gd.maize.lines}},
  \code{\link{gd.maize.aflp}}
}

\keyword{datasets}
