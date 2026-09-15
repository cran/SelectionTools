\name{plant}
\alias{plant}
\title{Create or process a plant object in the simulation framework}

\description{
Create or process a plant object in the simulation framework.
}

\usage{
plant(genotype, NoInd=1)
}

\arguments{
   \item{genotype}{Function argument used by this operation.}
   \item{NoInd}{Function argument used by this operation.}
}
\value{
A \code{data.frame} with columns \code{ind}, \code{chrom}, \code{hom}, \code{pos}, and \code{all}, representing the requested genome from the supplied homologue genotypes.
}
