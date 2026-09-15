\name{sel.target.define}
\alias{sel.target.define}
\title{Define a selection target}

\description{
Define a selection target.
}

\usage{
sel.target.define(criteria.df)
}

\arguments{
   \item{criteria.df}{Function argument used by this operation.}
}
\value{
A character vector containing one or two selection-target codes selected from \code{"H1"}, \code{"MIX"}, \code{"H2"}, and \code{"All"}. Returns \code{NULL} when the criterion cannot be translated into a supported target.
}
