\name{st.return.performance.data}
\alias{st.return.performance.data}
\title{Returns performance data from a SelectionTools data set}

\description{
Returns performance data from a SelectionTools data set.
}

\usage{
st.return.performance.data(out.filename = "y.vector", auxfiles = TRUE,
                           data.set = "default")
}

\arguments{
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a \code{data.frame} with character column \code{i}
(individual identifier) and numeric column \code{y} (performance value) when
the data can be returned; otherwise \code{NULL}. Numeric-looking individual
identifiers are retained as character strings.
}
