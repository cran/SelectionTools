\name{gd.allow.zero.frequencies}

\alias{gd.allow.zero.frequencies}

\title{gd.allow.zero.frequencies()}

\description{
  If for a marker no allele was observed, this is a missing value if
  codominant markers are used. However, for dominant markers where each
  band is regarded as a marker with only one allele, markers
  where no allele are observed are not missing data.
}

\usage{
gd.allow.zero.frequencies(allow=0)
}

\arguments{
   \item{allow}{Codominant marker data (SSR or RFLP) is analyzed like the datasets
  \code{gd.maize.lines} or \code{gd.maize.populations}. Dominant marker data (AFLP) is analyzed like the dataset
  \code{gd.maize.aflp}.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of setting whether zero allele frequencies are allowed in the genetic-diversity routines.
}
