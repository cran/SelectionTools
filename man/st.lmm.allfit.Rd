\name{st.lmm.allfit}
\alias{st.lmm.allfit}
\alias{summary.st.lmm.allfit}
\alias{print.st.lmm.allfit}
\title{Compares Gaussian mixed-model optimizers and starting values}
\description{
Refits a compiled Gaussian model with the available optimizers and specified
starting values. Every attempt uses the same observations, fixed and random
designs, covariance families, likelihood definition and parameter bounds.
Successful, nonconverged and failed attempts are retained separately.
}
\usage{
st.lmm.allfit(object, optimizers = NULL,
  starts = c("fitted", "default"), control = list(), keep.fits = TRUE)

\method{summary}{st.lmm.allfit}(object, ...)
\method{print}{st.lmm.allfit}(x, ...)
}
\arguments{
\item{object}{An ML or REML \code{st.lmm} fit; for \code{summary}, the
result of \code{st.lmm.allfit}. The original fit need not have converged,
but each selected starting vector must define a finite Gaussian likelihood.}
\item{optimizers}{Unique selections from \code{"AI"} and \code{"L-BFGS-B"}.
The default compares both for REML and uses \code{"L-BFGS-B"} for ML.
Models with explicit covariance constraints or a US trait-by-occasion
residual product use only \code{"L-BFGS-B"};
an explicit AI request for such a model is rejected.
Average information is available for REML. Its safeguarded implementation
may continue with L-BFGS-B; the requested and actual optimizer and full
fallback diagnostics are recorded.}
\item{starts}{Unique selections of \code{"fitted"} and \code{"default"},
or a uniquely named list of numerical coordinate vectors in the order
returned by \code{st.lmm.objective}. Named vectors must match that exact
coordinate order. The default starting vector is constructed from the
retained model with estimated covariance starts cleared. Every selected
optimizer receives exactly the same vector for a given start.}
\item{control}{Common fitting controls. This list must not set
\code{optimizer}, \code{parameters} or \code{evaluate.only}; the first two
are determined by the comparison arguments.}

\item{keep.fits}{Logical. Full fitted objects are retained in addition to estimates and
diagnostics when TRUE.}
\item{x}{A \code{st.lmm.allfit} object.}
\item{...}{Display controls for the result-table printer. Omitted controls use
\code{\link{st.options}}.}
}
\details{
Repeated fits run in the current R process. Numerical kernels use the
existing OpenMP and BLAS thread controls; \code{\link{st.set.num.threads}} describes these controls.
Intermediate state stays in C-owned memory.

This function compares the optimizers implemented by \code{st.lmm}; it
does not call external optimization packages. Repeated fits reuse the
retained model and do not reread mutable data or a marker-data registry.
Subset selection and random-effect tests are disabled in these repeated
fits, and fixed-effect inference uses the asymptotic setting so that the
comparison measures fitting rather than additional inference procedures.
The supplied model is not changed or automatically replaced.

Explicit fixed/shared covariance constraints are retained for all attempts,
including compact results with \code{keep.fits=FALSE}. Each start uses the
same physical parameter declarations. Compact estimates retain the constraint
inventory and independent covariance-parameter count.

The summary uses the converged attempt with the highest log likelihood as
its reference. It reports the initial and final deviance, deviance
difference, projected gradient, covariance-cone score, rank, differences
in fixed estimates, root mean squared difference of fitted values, and
absolute and relative Frobenius differences of physical covariance-block
matrices. Fixed-coefficient differences are omitted if the retained
estimable coefficient basis changes. If no attempt converged, reference
comparisons are missing. A nonconverged attempt with a better objective
can have a negative deviance difference and must be investigated rather
than treated as a valid reference fit.

Warnings, error messages, timings, coordinates, bounds, estimates and
execution identities remain in individual task records.
}
\value{
A \code{st.lmm.allfit} object containing the original fit, task inventory,
initial objective records, all attempt records, controls, execution
metadata. \code{summary} returns a data frame with
one row per attempt. Its \code{reference.id} attribute identifies the
successful reference attempt, or is missing when none succeeded.
}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lmm.refit}},
\code{\link{st.lmm.objective}}, \code{\link{st.lmm.bootstrap}},
\code{\link{st.lmm.profile}}}
