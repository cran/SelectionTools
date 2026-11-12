\name{st.select.phen}
\alias{st.select.phen}
\title{Selects individuals using phenotypic information}

\description{
Selects individuals using phenotypic information supplied as a data frame or
stored in a SelectionTools data set.
}

\usage{
st.select.phen(p = NULL, n = 0, t = 0, decreasing = TRUE, nsmall = 1,
               data.set = NULL)
}

\arguments{
   \item{p}{Optional data frame containing individual identifiers in a column
   named \code{i} or \code{ii} and numeric performance values in a column named
   \code{y} or \code{yy}. If \code{p} is \code{NULL}, performance data are
   obtained from \code{data.set}.}
   \item{n}{Maximum number of individuals to return. A value of 0 returns all
   individuals remaining after threshold selection.}
   \item{t}{Optional performance threshold. If greater than 0, values greater
   than \code{t} are retained when \code{decreasing=TRUE}; values smaller than
   \code{t} are retained when \code{decreasing=FALSE}.}
   \item{decreasing}{Logical. If \code{TRUE}, sort performance from largest to
   smallest; otherwise sort from smallest to largest.}
   \item{nsmall}{Minimum number of digits after the decimal point used in the
   generated \code{descr} column.}
   \item{data.set}{Optional name of a SelectionTools data set. Used when
   \code{p=NULL}; the performance values stored in that data set are selected.}
}
\value{
A \code{data.frame} containing the selected individuals, sorted and filtered as
requested, with an added character column \code{descr} combining the individual
identifier and performance value. The original identifier/performance column
names and column types are retained when \code{p} is supplied. When performance
data are obtained from \code{data.set}, individual identifiers are character
strings.
}
