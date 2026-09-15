\name{population.matrix.save}
\alias{population.matrix.save}
\title{Matrix save a population or population data}

\description{
Matrix save a population or population data.
}

\usage{
population.matrix.save(file, PopNames, missing.rm=TRUE, keep.IndID=TRUE)
}

\arguments{
   \item{file}{File name or file connection.}
   \item{PopNames}{Population name or population specification.}
   \item{missing.rm}{Function argument used by this operation.}
   \item{keep.IndID}{Function argument used by this operation.}
}
\value{
\code{NULL}, returned invisibly by the file-writing operation. The function is called for its side effect of writing matrix-format population data.
}
