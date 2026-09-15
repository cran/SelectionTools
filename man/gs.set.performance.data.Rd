\name{gs.set.performance.data}
\alias{gs.set.performance.data}
\title{Set performance data used by the genomic-selection routines}

\description{
Set performance data used by the genomic-selection routines.
}

\usage{
gs.set.performance.data(y, data.set = "default")
}

\arguments{
   \item{y}{Response vector.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
An invisible integer scalar containing the status value returned by the compiled performance-data setter.
}
