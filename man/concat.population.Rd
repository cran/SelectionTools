\name{concat.population}
\alias{concat.population}
\title{concat.population()}

\description{
Concatenates two populations and removes the second population.
}

\usage{
concat.population(NameP1, NameP2)
}

\arguments{
   \item{NameP1}{Recipient population}
   \item{NameP2}{Name of the population to be concatenated}
}
\value{
Returns the same invisible implementation-level \code{list} as \code{population.concat()}. The function is a compatibility wrapper and the meaningful result is the delegated population operation.
}

\details{
The individuals of \code{NameP2} are appended to \code{NameP1}. After successful concatenation, \code{NameP2} is removed. Complete individual records are moved or copied as required so chromosome data and any valid per-individual derived information remain associated with the correct individuals.
}

