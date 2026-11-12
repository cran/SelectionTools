\name{v-lmm-climate}
\alias{v-lmm-climate}
\alias{climate}
\docType{data}
\title{Balanced Climate and Variety Experiment}
\description{A balanced experiment with two varieties, three climates and four observations in each variety by climate combination.}
\usage{data("v-lmm-climate", package = "SelectionTools")}
\format{A \code{data.frame} named \code{climate} with 24 rows and 5 columns.
\describe{
\item{\code{variety}}{Variety identifier. Factor with levels a and b.}
\item{\code{climate}}{Climate identifier. Factor with levels I, II and III.}
\item{\code{cv}}{Variety by climate combination. Factor with six levels.}
\item{\code{rep}}{Replicate identifier. Factor with levels r1 to r4.}
\item{\code{yield}}{Yield response. Integer.}
}}
\details{The stored data contain no missing values.}
\value{Loading the dataset creates the \code{climate} data frame.}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lsm}}, \code{\link{v-lmm-data}}}
\keyword{datasets}
