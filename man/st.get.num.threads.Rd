\name{st.get.num.threads}
\alias{st.get.num.threads}
\title{Reports the computational thread configuration}
\description{
Returns the package thread budget, or report the package, OpenMP and BLAS
configuration without launching an OpenMP team.
}
\usage{
st.get.num.threads(complete = FALSE)
}
\arguments{
  \item{complete}{If \code{FALSE}, return the package thread budget as one
  integer. If \code{TRUE}, return the detailed configuration as a named list.}
}
\details{
The default result is the stored budget for SelectionTools OpenMP kernels;
it is not a measurement of an active team. It is 1 in a build without OpenMP.
Runtime and backend limits and operation size may reduce actual team sizes.
The startup hook calls \code{st.set.num.threads(2)}. An explicit setter call
caps the effective OpenMP budget at the runtime thread limit, disables
dynamic team shrinking, and applies the effective budget to supported BLAS
controls. A build without OpenMP limits an OpenMP OpenBLAS, or FlexiBLAS
with an unknown backend threading model, to one BLAS thread.

BLAS detection uses the library supplying R's linked BLAS operations.
If that library does not expose a supported control interface, the reported
BLAS count is \code{NA} and \code{blas.control} is \code{FALSE}. This includes
the serial reference BLAS supplied with R: the absence of a query interface
is not treated as evidence that an unidentified BLAS is serial.
}
\value{
With \code{complete = FALSE}, an integer scalar between 1 and 256.
With \code{complete = TRUE}, a named list containing:
\describe{
  \item{requested}{The initial budget or most recent explicit setting,
  after clamping to the range 1 to 256.}
  \item{package.threads}{The effective budget used by package OpenMP
  kernels, capped by the runtime thread limit, or 1 without OpenMP.}
  \item{explicit}{Whether a setter has been called in this process.}
  \item{openmp.available}{Whether SelectionTools was built with OpenMP.}
  \item{openmp.max.threads}{The OpenMP maximum for the calling thread,
  or 1 without OpenMP.}
  \item{openmp.dynamic}{Whether the OpenMP runtime may adjust team sizes,
  or \code{FALSE} without OpenMP.}
  \item{openmp.thread.limit}{The runtime's limit on OpenMP threads,
  or 1 without OpenMP.}
  \item{blas.backend}{The identified interface: \code{"OpenBLAS"},
  \code{"FlexiBLAS"}, \code{"MKL"}, or \code{"unavailable"}.}
  \item{blas.library}{The path of the library supplying the linked BLAS,
  when available.}
  \item{blas.threads}{The configured thread count reported by the backend,
  or \code{NA} when it cannot be queried. It is not a measurement of the
  threads used by a particular BLAS operation.}
  \item{blas.control}{Whether the supported backend's thread setter and
  getter are available.}
  \item{blas.parallel}{The OpenBLAS threading mode: 0 for a serial build,
  1 for platform threads, or 2 for OpenMP; \code{NA} when unavailable.}
  \item{legacy.outer.threads}{The thread count for older genomic-selection
  outer loops that call BLAS. These loops use 1 to avoid unsafe concurrent
  BLAS calls; BLAS itself retains its configured budget.}
}
}
\seealso{
\code{\link{st.set.num.threads}}
}
