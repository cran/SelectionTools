\name{st.lmm.profile}
\alias{st.lmm.profile}
\alias{profile.st.lmm}
\alias{confint.st.lmm.profile}
\alias{print.st.lmm.profile}
\alias{plot.st.lmm.profile}
\title{ML Likelihood Profiles for Gaussian Mixed Models}
\description{
Constructs likelihood profiles for fixed coefficients and supported random
and residual standard deviations and correlations.
All other estimable parameters are re-estimated at each specified target.
Confidence limits and signed-root likelihood-ratio plots use the resulting
ML profile, including an ML refit when the original model used REML.
}
\usage{
st.lmm.profile(object, parm = NULL, level = .95, control = list())

\method{profile}{st.lmm}(fitted, ...)

\method{confint}{st.lmm.profile}(object, parm,
    level = object$level, ...)

\method{print}{st.lmm.profile}(x, ...)

\method{plot}{st.lmm.profile}(x, which = NULL, ...)
}
\arguments{
\item{object,fitted}{An \code{st.lmm} model fitted by ML or REML for
\code{st.lmm.profile} and \code{profile}; an \code{st.lmm.profile} result
for its confidence-interval method. Method-of-moments fits have no likelihood
profile and are rejected.}
\item{parm}{Parameter names or integer indices. For \code{st.lmm.profile},
\code{NULL} selects all supported parameters. \code{parm = character()} is available
to obtain the parameter inventory without tracing any curves; its
\code{parameters} component gives names, estimates, physical bounds and
identifiability. For \code{confint} applied to a profile result, a missing
\code{parm} selects all previously profiled parameters, and indices refer to
that selected set. Duplicate and unknown parameters are errors.}
\item{level}{A confidence level strictly between zero and one. Requesting a
different level from \code{confint} on a saved profile retraces the requested
parameters at that level.}
\item{control}{Named list of profile controls. The defaults are
\code{maxit = 1000}, \code{maxsteps = 24}, \code{maxroot = 60},
\code{maxrestarts = 8}, \code{score.tol = 1e-7},
\code{cone.tol = 1e-7}, \code{lr.tol = 1e-5},
\code{constraint.tol = 1e-10}, and \code{extent = 1e6}.
\code{maxit} limits the initial nuisance-optimization iterations per start
and restart; a subsequent Newton refinement can take up to 35 steps.
\code{maxit = 0} evaluates the supplied starts without optimization.
\code{maxsteps} limits outward bracketing steps per side;
\code{maxroot} limits endpoint refinement; \code{maxrestarts} limits
restarts from feasible covariance-boundary directions.
\code{score.tol} and \code{cone.tol} bound the nuisance projected score
and covariance-cone violation. \code{lr.tol} bounds the error in the endpoint
likelihood-ratio statistic. \code{constraint.tol} bounds the relative error
in the fixed physical target. These four tolerances may be tightened, but
not relaxed beyond their defaults. \code{extent} limits the absolute search
distance to its value times \code{max(1, abs(estimate))}.
Unknown controls are errors.}

\item{x}{An \code{st.lmm.profile} result.}
\item{which}{Names or indices of profile curves to plot. \code{NULL} plots
all selected curves in their stored order.}
\item{...}{For printing, display controls whose omitted values use
\code{\link{st.options}}. Arguments passed from \code{profile} to \code{st.lmm.profile},
or named graphical arguments for \code{plot}. The profile \code{print} and
\code{confint} methods accept the generic's remaining arguments but do not
use them.}
}
\details{
Likelihood profiles for separable covariance products, including the
established AR1 by AR1 family, and US trait-by-occasion products are not implemented. Such requests raise an
explicit error; no unconstrained covariance model is substituted.

Repeated fits run in the current R process. Numerical kernels use the
existing OpenMP and BLAS thread controls; \code{\link{st.set.num.threads}} describes these controls.
R callbacks execute on the main R thread. Intermediate state stays in memory.

Parameter names begin with \code{beta:}, \code{sd:}, or \code{cor:}.
The remaining labels identify the covariance block and its coefficient or
coefficient pair. Standard deviations are marginal standard deviations in
the original coefficient units. They are not Cholesky coordinates or
standard deviations divided by the residual scale.

The inventory includes estimable fixed coefficients, residual standard
deviations, independent or unstructured random-effect standard deviations,
and every pairwise correlation of an unstructured random block. Homogeneous
random or residual covariance structures contribute one common standard
deviation; heterogeneous structures contribute one per coefficient or
occasion. Compound-symmetric and AR1 random or residual covariance
structures contribute their common correlation
parameter, subject to the corresponding dimension-dependent or unit bounds.
Fixed covariance blocks do not contribute estimated parameters. Pairwise
correlations of an unstructured residual block are optimized as nuisance
parameters but are not available as individual profile targets.

The calculation first obtains a converged ML baseline with the retained
response, observation mask, fixed and random designs, precision weights,
offsets, and covariance specification. This uses the retained resolved model
also for SelectionTools marker and hybrid data sets. The original fitted
object is unchanged. The result records both \code{method = "ML"} and
\code{original.method}; REML-based small-sample degrees of freedom are not
used to calibrate the profile.
\code{confint(fit, method = "profile", parm = ...)} provides a direct interval
interface. The existing \code{confint(fit)} default remains the Wald interval
for fixed coefficients.

Every curve uses the same ML baseline. Numerical records retain warnings
and failed nuisance-optimization starts. A materially better ML optimum
invalidates the requested intervals; numerical details remain available
for diagnosis.

For a parameter value \eqn{t}, the likelihood-ratio statistic is twice the
ML log-likelihood difference between the baseline and the constrained fit.
The confidence limit solves the equation equating this statistic to the
\code{level} quantile of a chi-squared distribution with one degree of
freedom. Its signed square root is negative below the estimate and positive
above it. This likelihood-ratio calibration is approximate in small samples
and at variance or correlation boundaries.

Each constrained fit compares the preceding successful point, the fitted
baseline, and a cold start. Positive-semidefinite covariance constraints and
the selected target are preserved throughout nuisance optimization.
Convergence requires the nuisance score and feasible covariance-direction
checks, as well as the target constraint. Tiny free structured standard
deviations may be projected to zero as a numerical restart only when the
objective change is within a recorded floating-point roundoff bound and
the target is preserved. All convergence checks are then repeated. The
projections and every start are retained in the result. A critical-cone
calculation involving more than 200 nuisance covariance coordinates, or an
unresolved degenerate cone within its work limit, yields an incomplete
point rather than an accepted interval.
The same conservative result applies to a partial-zero heterogeneous AR1
residual block when its correlation is unidentified and boundary curvature
over the equivalent correlations remains uncertified.

A pairwise correlation of an unstructured random block is unidentified
when either associated standard deviation is zero. The shared correlation
of a compound-symmetric or AR1 block requires at least one informative pair
with positive standard deviations. For residual blocks, those positions
must also be observed together within a group. If informative AR1 pairs
have only even lags, the sign of the correlation is unidentified. Aliased
fixed coefficients are also marked unidentified. These requests remain
visible in the output and do not yield ordinary confidence limits.

The endpoint statuses distinguish \code{complete} (the likelihood-ratio
threshold was reached), \code{boundary} (the physical bound was reached
within the accepted set), \code{flat} (approximately flat over the examined
range), \code{open} (the search limit was reached without a bracket),
\code{unidentified}, and \code{failed}. An open search does not by itself
establish an infinite confidence limit. \code{confint} returns \code{NA}
for an endpoint that is not \code{complete} or \code{boundary}, and attaches
a matrix of endpoint statuses. A constrained start that achieves a
materially better ML likelihood than the baseline triggers a warning and
invalidates all requested intervals; all start records remain available for
diagnosis. A lower objective from an unconverged start cannot be silently
discarded in favour of a worse converged start.

Profile plotting uses the ordinary graphics layout and does not change
\code{par(mfrow)} or force a multi-panel arrangement. Each requested curve
is plotted with horizontal reference lines at the confidence threshold.
}
\value{
\code{st.lmm.profile} and \code{profile} return a serializable object of
class \code{st.lmm.profile}. Its principal components are:
\itemize{
\item \code{parameters}: the complete physical parameter inventory.
\item \code{selected}: names of the requested parameters.
\item \code{intervals}: estimates, endpoints, endpoint statuses,
likelihood-ratio errors, searched values and point identifiers.
\item \code{profiles}: a concise point table and detailed start records for
each selected parameter. The detailed records contain physical covariance
matrices, fixed coefficients, nuisance scores, covariance-cone diagnostics,
constraint errors, numerical coordinates and their exact constrained map,
and any accepted boundary projections.
\item \code{baseline}, \code{baseline.deviance}, \code{method}, and
\code{original.method}: the ML reference fit and its provenance.
\item \code{control}, \code{level}, \code{cutoff}, \code{complete}, and
\code{diagnostic}: the effective settings and overall completion status.
\item \code{object}: the original fitted model, retained for retracing.
}
The confidence-interval method returns a two-column numeric matrix with
\code{status}, \code{method} and \code{level} attributes. The print and plot
methods return their profile object invisibly.
}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lmm.refit}},
\code{\link{st.lmm.objective}}, \code{\link[stats]{confint}}}
