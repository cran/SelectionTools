\name{v-lmm-tomato}
\alias{v-lmm-tomato}
\alias{tomato}
\docType{data}
\title{Tomato Plant Heights for Linear Hypotheses}
\description{Seven plant-height observations in three treatment groups, used for a cell-means model and explicit linear hypotheses.}
\usage{data("v-lmm-tomato", package = "SelectionTools")}
\format{A \code{data.frame} named \code{tomato} with 7 rows and 2 columns.
\describe{
\item{\code{y}}{Plant height. Numeric (double).}
\item{\code{x.v}}{Treatment identifier. Factor with levels 1, 2 and 3; group sizes are 3, 2 and 2.}
}}
\details{Plant-height measurement units are not recorded.
The stored data contain no missing values.}
\value{Loading the dataset creates the \code{tomato} data frame.}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lsm}}, \code{\link{v-lmm-data}}}
\keyword{datasets}
