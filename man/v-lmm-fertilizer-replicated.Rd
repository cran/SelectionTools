\name{v-lmm-fertilizer-replicated}
\alias{v-lmm-fertilizer-replicated}
\alias{fertilizer.replicated}
\docType{data}
\title{Replicated Fertilizer Doses for a Regression Lack-of-Fit Test}
\description{Eleven yield observations at four fertilizer doses, used to
partition regression residual variation into lack of fit and pure error.}
\usage{data("v-lmm-fertilizer-replicated", package = "SelectionTools")}
\format{A \code{data.frame} named \code{fertilizer.replicated} with 11 rows
and 2 columns.
\describe{
\item{\code{dueng}}{Numeric fertilizer dose, with values 1, 3, 5 and 7.
The first three doses have three observations each; dose 7 has two.}
\item{\code{ertrag}}{Numeric crop yield.}
}}
\details{Measurement units are not recorded. There are no missing values.
Replication at the four dose levels provides 7 degrees of freedom for
pure error. A linear mean function has 2 degrees of freedom for lack of fit
relative to unrestricted dose means.}
\value{Loading the dataset creates the \code{fertilizer.replicated}
data frame.}
\keyword{datasets}
