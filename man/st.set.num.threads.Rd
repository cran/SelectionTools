\name{st.set.num.threads}
\alias{st.set.num.threads}
\title{Sets the computational thread budget}
\description{
Sets the thread budget for SelectionTools OpenMP calculations and, when its
runtime interface is available, for the BLAS library used by R.
}
\usage{
st.set.num.threads(active)
}
\arguments{
  \item{active}{One finite integer. Values below 1 are set to 1 and values
  above 256 are set to 256.}
}
\details{
The explicit setting is stored by SelectionTools and supplied to its parallel
numerical kernels. It also sets the calling thread's OpenMP maximum and the
configured thread count of a supported BLAS. OpenBLAS, FlexiBLAS and MKL
interfaces are used only when their control routines can be identified in the
library supplying R's linked BLAS operations. No BLAS library or R package is
loaded to provide these interfaces. Other BLAS implementations remain usable,
but their thread counts cannot be queried or changed by this function.

The package startup hook calls this function with a budget of 2. A subsequent call
after loading the package selects another budget. Older routines
allocate per-thread work for at most 256 threads, so the requested value is
clamped to that range. With OpenMP, the effective package and BLAS budgets
are also capped by the runtime thread limit. An explicit setting disables
OpenMP dynamic team shrinking to keep the package and OpenMP BLAS budgets
consistent. The setting applies to the current R process; other packages
using the same OpenMP runtime or BLAS may also be affected.

Backend limits and the size of an operation can reduce the number of threads
actually used. Builds without OpenMP perform package calculations serially.
A detected pthread OpenBLAS can still use the requested BLAS budget. A
build without OpenMP limits an OpenMP OpenBLAS, or FlexiBLAS with an unknown
backend threading model, to one BLAS thread because it cannot synchronize
the OpenMP runtime. This avoids waiting for BLAS workers that cannot start.
The reference BLAS supplied with R is serial and has no thread-control
interface.

The parallel kernels used by \code{st.lmm} do not call BLAS or the R API
from their worker threads. The older genomic-selection cross-validation and
permutation-test loops call BLAS and therefore execute their outer iterations
serially. This avoids concurrent BLAS calls from OpenMP teams with backends
for which that combination is unsafe. BLAS operations in these iterations
retain the configured BLAS budget, so the operations themselves can still
use multiple threads. This conservative policy can reduce parallel speedup
of those older outer loops.

\code{st.get.num.threads(complete = TRUE)} reports the stored budget,
OpenMP runtime settings and available BLAS information. Environment-variable
values alone do not verify the current runtime settings.
}
\value{
An invisible list containing the clamped integer setting returned by the
compiled routine. The function is called for its effect on the
thread configuration.
}
\seealso{
\code{\link{st.get.num.threads}}, \code{\link{gs.set.num.threads}}
}
