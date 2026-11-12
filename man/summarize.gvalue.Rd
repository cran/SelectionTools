\name{summarize.gvalue}
\alias{summarize.gvalue}
\title{summarize.gvalue()}

\description{
Summarizes the genetic values of the individuals of a population with
respect to defined effects
}

\usage{
summarize.gvalue(pops, effectfile=NULL)
}

\arguments{
   \item{pops}{Names separated by blanks of the populations to be
   summarized.)}
   \item{effectfile}{Name of the effect}
}
\value{
A numeric \code{matrix} with one row per population and columns \code{Obs}, \code{Mean}, \code{SDev}, \code{Min}, \code{Q10}, \code{Med}, and \code{Max}, summarizing genetic values.
}

\details{
If no effects are specified the weighted sum of all effects are used,
otherwise only the single specified effects is used without using weights.
}

\note{Further details appear in the \code{lib00.example1} for an application of the command.}
