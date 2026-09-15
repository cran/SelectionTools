\name{gd.list.irregular}

\alias{gd.list.irregular}

\title{gd.list.irregular()}

\description{
  The function \emph{gd.list.irregular} lists individuals for which at one
  marker more than two alleles were observed.
}

\usage{
gd.list.irregular(dta1)
}

\arguments{
   \item{dta1}{Name of the dataset containing the observed marker
   genotypes. The format of the dataset is illustrated by the datasets
 gd.dummy.populations and gd.maize.populations.}
}

\value{
A \code{data.frame} with columns \code{ind} and \code{marker}, identifying individual/marker combinations for which more than two alleles were observed.
}
