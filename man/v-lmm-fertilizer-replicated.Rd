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
For a lack-of-fit test, include \code{dueng} followed by \code{factor(dueng)}
in the fixed formula and request sequential tests with
\code{method="moments"} and \code{type=1}.
The \code{factor(dueng)} test returned by \code{anova} compares the linear
mean function with unrestricted dose means, using the pure-error residual.
Pure error has 7 df and lack of fit has 2 df.}
\value{Loading the dataset creates the \code{fertilizer.replicated}
data frame.}
\keyword{datasets}
