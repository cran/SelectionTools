\name{gs.set.performance.data}
\alias{gs.set.performance.data}
\title{Sets performance data used by the genomic-selection routines}

\description{
Sets performance data used by the genomic-selection routines.
}

\usage{
gs.set.performance.data(y, data.set = "default")
}

\arguments{
   \item{y}{Numeric response vector in the order of the individuals in the
   data set. Its length must equal the number of individuals. Missing values
   (\code{NA} or \code{NaN}) can mark individuals without observations;
   infinite values are rejected.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
An invisible integer scalar containing the status value returned by the compiled performance-data setter.
}
\details{
All individuals and marker data are retained. The stored vector can be
retrieved with \code{gs.get.y}. The mixed-model interface \code{st.lmm}
supports missing responses. Routines that require complete response vectors
must be given only the observed values and corresponding inputs.
}
