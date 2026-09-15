\name{st.return.performance.data}
\alias{st.return.performance.data}
\title{Return performance data from a SelectionTools data set}

\description{
Return performance data from a SelectionTools data set.
}

\usage{
st.return.performance.data(out.filename = "y.vector", auxfiles = TRUE,
                           data.set ="default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} with columns \code{i} (individual identifier) and \code{y} (performance value) when the data can be returned; otherwise \code{NULL}.
}
