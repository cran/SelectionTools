\name{gd.correct.missing}

\alias{gd.correct.missing}

\title{gd.correct.missing()}

\description{
  The function \emph{gd.correct.missing} corrects missing marker data
  points in a data set. Two rules are applied: (1) If at least one
  allele was observed at a marker, missing values for this marker at
  other alleles are changed to zeros. (2) If at a marker no allele was
  observed but for all alleles only zeros are in the data set,
  the zeros are changed to missing data.
}

\usage{
gd.correct.missing(dta)
}

\arguments{
   \item{dta}{Name of the dataset contaning the observed marker
   genotypes. The format of the dataset is illustrated by the datasets
 gd.dummy.populations and gd.maize.populations.}
}

\value{
A \code{data.frame} containing the corrected data, with the same row and column structure as the supplied data object.
}
