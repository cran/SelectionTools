\name{select.n.best.segments}
\alias{select.n.best.segments}
\title{Selects the best genome segments according to the requested criterion}

\description{
Selects the best genome segments according to the requested criterion.
}

\usage{
select.n.best.segments(newPop, oldPop, n=1, allele = 1, decreasing = FALSE)
}

\arguments{
   \item{newPop}{Function argument used by this operation.}
   \item{oldPop}{Function argument used by this operation.}
   \item{n}{Number of items or repetitions.}
   \item{allele}{Allele specification.}
   \item{decreasing}{Function argument used by this operation.}
}
\value{
On successful selection, returns the invisible implementation-level \code{list} produced by \code{population.divide()}; if the requested number cannot be selected, \code{NULL} is returned invisibly. The main effect is creation of \code{newPop} from individuals ranked by the number of genome segments.
}

\details{
Each derived cache is complete for all individuals or absent throughout a
population. A mixed source cache causes all marker-genotype, genetic-value and
phenotype caches to be cleared silently in the source and result.
}
