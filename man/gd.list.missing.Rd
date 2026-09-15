\name{gd.list.missing}

\alias{gd.list.missing}

\title{gd.list.missing()}

\description{
  The function \emph{gd.list.missing} lists the missing marker data points
  in a data set containing molecular marker data.
}

\usage{
gd.list.missing(dta)
}

\arguments{
   \item{dta}{Name of the dataset containing the observed marker
   genotypes. The format of the dataset is illustrated by the datasets
 gd.dummy.populations and gd.maize.populations.}
}

\value{
A \code{data.frame} with columns \code{ind} and \code{marker}, identifying individual/marker combinations with missing observations.
}
