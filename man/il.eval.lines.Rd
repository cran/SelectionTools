\name{il.eval.lines}
\alias{il.eval.lines}
\title{Constructs or evaluates introgression-line libraries}

\description{
Constructs or evaluates introgression-line libraries.
}

\usage{
il.eval.lines(allele = 1, chrom = 0, lower = 0, upper = 0, exclude = 0,
              data.set = "default")
}

\arguments{
   \item{allele}{Allele specification.}
   \item{chrom}{Chromosome identifier.}
   \item{lower}{Function argument used by this operation.}
   \item{upper}{Function argument used by this operation.}
   \item{exclude}{Function argument used by this operation.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
A \code{data.frame} with one row per evaluated line and columns \code{name}, \code{nseg}, \code{lseg}, \code{abs}, \code{tot}, and \code{rel}, summarizing introgressed segments and their lengths/proportions. Returns \code{NULL} on validation failure.
}
