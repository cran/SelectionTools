\name{generate.population}
\alias{generate.population}
\title{generate.population()}

\description{
Generates one or more simulation populations from the marker-incidence matrix representation used by the genetic-distance routines.
}

\usage{
generate.population(dta, backcross=FALSE)
}

\arguments{
\item{dta}{Data frame or matrix in the marker-incidence representation expected by the genetic-distance population interface.}
   \item{backcross}{Logical flag selecting backcross-oriented generation in the backend.}
}
\value{
An invisible implementation-level \code{list} returned by either the population-generation routine or the information routine used to report invalid allele encoding. The meaningful result is the side effect of generating the requested population when the input is valid.
}

\details{
The function interprets the row and column descriptors of \code{dta} to determine marker names, allele codes, population names, and individual membership. Allele coding must be numeric integer coding compatible with the simulation routines.

The compiled interface derives the population and marker dimensions and supplies the matrix together with allele and missing-value information to the population generator. The resulting objects are simulation populations represented by chromosome data, not SelectionTools marker-data sets.

The \code{backcross} flag selects the corresponding population-generation mode in the simulation backend. Existing genome and linkage-map definitions determine how generated marker information is represented in simulated chromosomes.
}
