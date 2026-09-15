\name{ph.linkage.map.create}
\alias{ph.linkage.map.create}
\title{Create a linkage map for phenotype or simulation use}

\description{
Create a linkage map for phenotype or simulation use.
}

\usage{
ph.linkage.map.create(map, disperse=FALSE, file.disperse=NA,
                      disperse.factor=100)
}

\arguments{
   \item{map}{Linkage map or map object.}
   \item{disperse}{Function argument used by this operation.}
   \item{file.disperse}{Function argument used by this operation.}
   \item{disperse.factor}{Function argument used by this operation.}
}
\value{
When \code{disperse = TRUE}, returns the modified map object after positions have been dispersed. With the default \code{disperse = FALSE}, returns \code{NULL}; in both cases the main side effect is installing the map through the compiled map-definition routine.
}
