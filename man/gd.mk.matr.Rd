\name{gd.mk.matr}
\alias{gd.mk.matr}
\title{Construct a matrix from genetic-distance results}

\description{
Construct a matrix from genetic-distance results.
}

\usage{
gd.mk.matr(dta, distance)
}

\arguments{
   \item{dta}{Input data object.}
   \item{distance}{Distance specification or matrix.}
}
\value{
A numeric \code{matrix} representation of the supplied pairwise distance results, with population labels taken from the data object.
}
