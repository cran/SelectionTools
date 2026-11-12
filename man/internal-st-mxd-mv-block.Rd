\name{.st_mxd_mv_block}
\alias{.st_mxd_mv_block}
\title{Constructs a Covariance-Block Specification}
\description{
Constructs a covariance-block description for the matrix-based Gaussian
mixed-model engine.
}
\usage{
.st_mxd_mv_block(group = NULL, basis = NULL, A = NULL,
                 covariance = "diag", start = NULL, kernel = NULL, hom = NULL)
}
\arguments{
\item{group}{Positive integer group indices, a factor or character group
labels, with one entry per response cell for a grouped block.}
\item{basis}{A numeric coefficient-basis matrix with one row per response
cell. A vector is treated as a single column. If \code{NULL}, a column of
ones is used for a grouped block.}
\item{A}{An optional relationship matrix among group levels. Its rows and
columns must cover the groups; supplied labels are used for alignment.}
\item{covariance}{The coefficient covariance structure: \code{"diag"},
\code{"unstructured"}, \code{"fixed"}, \code{"homdiag"}, \code{"cs"} or
\code{"ar1"}. The last three structures are supported for univariate
random-coefficient blocks.}
\item{start}{Optional starting covariance values, or the supplied covariance
for a fixed block.}
\item{kernel}{An optional covariance kernel over response cells. It cannot
be combined with \code{group}, \code{basis} or \code{A}.}
\item{hom}{Optional logical value specifying a common coefficient variance
for \code{"cs"} or \code{"ar1"}. With \code{NULL}, compound symmetry
uses heterogeneous variances and first-order autoregression uses a common
variance. The \code{"homdiag"} structure always uses a common variance.
It is used only by these three covariance structures.}
}
\value{
A list describing a grouped covariance block or a response-cell kernel
block. Grouped specifications contain group indices, the coefficient
basis, relationship information, covariance settings and available level
labels. The fitting engine performs the numerical covariance validation.
}
\keyword{internal}
