\name{gd.allele.frequencies}

\alias{gd.allele.frequencies}

\title{gd.allele.frequencies()}

\description{
 Calculates the allele frequencies of the populations in a data set
 consisting of codominant molecular marker data.
}

\usage{
gd.allele.frequencies(dta)
}

\arguments{
   \item{dta}{Name of the dataset contaning the observed marker
   genotypes. The format of the dataset is illustrated by the datasets
   gd.dummy.populations and gd.maize.populations.}
}

\value{
A \code{data.frame} of allele frequencies, with marker/allele combinations in rows and populations in columns.
}
