\name{v-lmm-iris}
\alias{v-lmm-iris}
\alias{iris.manova}
\docType{data}
\title{Iris Measurements for MANOVA}
\description{Four measurements on 150 Iris plants from three species, used for multivariate analysis of variance.}
\usage{data("v-lmm-iris", package = "SelectionTools")}
\format{A \code{data.frame} named \code{iris.manova} with 150 rows and 5 columns.
\describe{
\item{\code{SepalLength}}{Sepal length in centimetres. Numeric (double).}
\item{\code{SepalWidth}}{Sepal width in centimetres. Numeric (double).}
\item{\code{PetalLength}}{Petal length in centimetres. Numeric (double).}
\item{\code{PetalWidth}}{Petal width in centimetres. Numeric (double).}
\item{\code{Species}}{Species identifier. Factor with levels setosa, versicolor and virginica.}
}}
\details{The measurements agree with the standard Iris dataset.
The stored data contain no missing values.}
\value{Loading the dataset creates the \code{iris.manova} data frame.}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lsm}}, \code{\link{v-lmm-data}}}
\keyword{datasets}
