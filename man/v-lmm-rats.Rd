\name{v-lmm-rats}
\alias{v-lmm-rats}
\alias{rats}
\docType{data}
\title{Rats for Partial and Sequential Regression Tests}
\description{Initial weight, feed intake and final weight for 35 rats, used to distinguish partial from sequential tests in multiple regression.}
\usage{data("v-lmm-rats", package = "SelectionTools")}
\format{A \code{data.frame} named \code{rats} with 35 rows and 3 columns.
\describe{
\item{\code{x.1}}{Initial rat weight. Numeric (double).}
\item{\code{x.2}}{Feed intake. Integer.}
\item{\code{y}}{Final rat weight. Numeric (double).}
}}
\details{Measurement units are not recorded.
The stored data contain no missing values.}
\value{Loading the dataset creates the \code{rats} data frame.}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lsm}}, \code{\link{v-lmm-data}}}
\keyword{datasets}
