\name{population.transfer}
\alias{population.transfer}
\title{Transfer a population or population data}

\description{
Transfers the individuals of one simulation population into another, either retaining or deleting the source population.
}

\usage{
population.transfer(NameP1, NameP2, population.NameP2.delete=TRUE)
}

\arguments{
\item{NameP1}{Character string naming the recipient population.}
   \item{NameP2}{Character string naming the source population.}
   \item{population.NameP2.delete}{Logical. If true, remove the source after transfer by concatenation; if false, retain the source and append a copy.}
}
\value{
The invisible implementation-level \code{list} returned by \code{population.concat()} or \code{population.append()}, depending on \code{population.NameP2.delete}. The main effect is transfer of population data.
}

\details{
If \code{population.NameP2.delete=TRUE}, the operation uses population concatenation: the individuals of \code{NameP2} are added to \code{NameP1} and \code{NameP2} is removed. If false, population append is used: the same individuals are copied to \code{NameP1} while \code{NameP2} remains available. Derived genotype and genetic-value state follows the consistency rules of the underlying append or concatenate operation.
}
