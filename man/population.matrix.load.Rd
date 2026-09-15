\name{population.matrix.load}
\alias{population.matrix.load}
\title{Matrix load a population or population data}

\description{
Matrix load a population or population data.
}

\usage{
population.matrix.load(file, backcross=FALSE, keep.IndID=TRUE)
}

\arguments{
   \item{file}{File name or file connection.}
   \item{backcross}{Function argument used by this operation.}
   \item{keep.IndID}{Function argument used by this operation.}
}
\value{
\code{NULL}. The function is called for its side effect of loading matrix-format population data and, when requested, restoring individual identifiers.
}
