\name{single.cross}
\alias{single.cross}
\title{Create or evaluate a single cross}

\description{
Create or evaluate a single cross.
}

\usage{
single.cross(PgName, PopName, classSize=1)
}

\arguments{
   \item{PgName}{Function argument used by this operation.}
   \item{PopName}{Name of the population.}
   \item{classSize}{Function argument used by this operation.}
}
\value{
An invisible implementation-level \code{list} returned by the final population-removal operation. The meaningful result is the side effect of creating the requested single-cross population.
}
