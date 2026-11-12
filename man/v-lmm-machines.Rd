\name{v-lmm-machines}
\alias{v-lmm-machines}
\alias{machines}
\docType{data}
\title{Machine and Operator Experiment}
\description{Performance scores from two machine types evaluated twice by
each of three operators. The crossed design illustrates inference spaces
for least-squares means in linear mixed models.}
\usage{data("v-lmm-machines", package = "SelectionTools")}
\format{A \code{data.frame} named \code{machines} with 12 rows and 3 columns.
\describe{
\item{\code{machine}}{Machine type. Factor with levels \code{1} and \code{2}.}
\item{\code{operator}}{Operator identifier. Factor with levels \code{1}, \code{2} and \code{3}.}
\item{\code{score}}{Observed performance score. Numeric (double).}
}}
\details{Each machine and operator combination has two observations.
The stored data contain no missing values.}
\value{Loading the dataset creates the \code{machines} data frame.}
\source{Sanders, W. L. (1991). The use of the SAS system in statistical training.
Proceedings of the Sixteenth Annual SAS Users Group International Conference.
\url{https://support.sas.com/resources/papers/proceedings-archive/SUGI91/Sugi-91-283\%20Sanders.pdf}}
\references{McLean, R. A., Sanders, W. L. and Stroup, W. W. (1991).
A unified approach to mixed linear models. The American Statistician, 45, 54--64.
\doi{10.1080/00031305.1991.10475767}.}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lsm}}, \code{\link{v-lmm-data}}}
\keyword{datasets}
