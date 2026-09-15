\name{st.select.phen}
\alias{st.select.phen}
\title{Select individuals using phenotypic information}

\description{
Select individuals using phenotypic information.
}

\usage{
st.select.phen(p, n=0, t=0, decreasing=TRUE, nsmall=1)
}

\arguments{
   \item{p}{Function argument used by this operation.}
   \item{n}{Number of items or repetitions.}
   \item{t}{Generation or method-specific time parameter.}
   \item{decreasing}{Function argument used by this operation.}
   \item{nsmall}{Function argument used by this operation.}
}
\value{
A \code{data.frame} containing the selected rows of \code{p}, sorted and filtered as requested, with an added character column \code{descr} combining individual identifier and phenotypic value.
}
