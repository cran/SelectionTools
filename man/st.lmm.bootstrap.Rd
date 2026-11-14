\name{st.lmm.bootstrap}
\alias{st.lmm.bootstrap}
\alias{summary.st.lmm.bootstrap}
\alias{print.st.lmm.bootstrap}
\alias{confint.st.lmm.bootstrap}
\title{Bootstrap Statistics from Gaussian Linear Mixed Models}
\description{
Generates responses from an ML or REML fit, refits its retained model design,
and evaluates a named numeric statistic. Replicate records include the
generated response, random-number stream, convergence diagnostics, warnings
and failures.
}
\usage{
st.lmm.bootstrap(object, statistic = stats::coef, nsim = 999L,
  seed = NULL, type = c("parametric", "residual"), re.form = NA,
  statistic.args = list(), responses = NULL, control = list(), series = NULL)

\method{summary}{st.lmm.bootstrap}(object, ...)
\method{print}{st.lmm.bootstrap}(x, ...)
\method{confint}{st.lmm.bootstrap}(object, parm, level = .95,
  type = c("percentile", "basic", "normal"), allow.failed = FALSE, ...)
}
\arguments{
\item{object}{For bootstrap generation, a converged \code{st.lmm} ML or
REML fit. For summary and intervals, an \code{st.lmm.bootstrap} result.}
\item{statistic}{A function whose first argument is a refitted model. It
must return a nonempty finite numeric vector with unique, nonempty names
in the same order for every replicate. The default returns fixed-effect
coefficients. A rank-deficient model requires a statistic selecting
estimable coefficients or combinations instead of returning aliased entries.}
\item{nsim}{Positive integer number of replicates.}
\item{seed}{An integer master seed, or \code{NULL} to draw one from R's
current random-number state. A supplied seed restores the caller's previous
RNG kind and state, including the absence of \code{.Random.seed}. With
\code{NULL}, drawing the recorded master seed advances the caller's state.}
\item{type}{For generation, \code{"parametric"} samples Gaussian random
effects and residuals; \code{"residual"} resamples centered and scaled
conditional residuals in the limited setting described below. For
\code{confint}, the interval construction.}
\item{re.form}{Random coefficients held at their fitted values during a
parametric bootstrap; \code{\link{simulate.st.lmm}} describes the generation step. The default
\code{NA} regenerates all random effects. \code{NULL} holds all of them.
Selectors may hold complete blocks or fitted coefficients within a block.
Residual resampling always holds every fitted random coefficient and
requires the default \code{re.form=NA}; it has no selectable conditioning
argument.}
\item{statistic.args}{A list of additional named arguments passed to
\code{statistic} after the fitted model.}
\item{responses}{Optional externally generated responses. A list of
response matrices has one matrix per simulation, with one row per retained
observation and one column per trait. A univariate model also accepts a list
of vectors or a matrix with one simulation per column. The number of
responses must equal \code{nsim}. Each response uses the retained layout of
\code{model.response(model.frame(object))}. Original missing response cells
remain missing; an observed cell cannot become nonfinite.}
\item{control}{Named fitting controls for repeated fits. The ordinary route
passes them to \code{\link{st.lmm.refit}}. The joint series route applies
the same controls to the compiled first- and second-stage refits; physical
\code{control$parameters} overrides are not supported there.}
\item{series}{NULL for the ordinary bootstrap, which holds any inherited
first-stage covariance fixed. Refitting both stages uses the original
earlier fits supplied as \code{list("env", fit1, fit2, ...)} in their original order.
The object must be a series fit with its complete original target cells.
Further details appear in \sQuote{Joint bootstrap of a series}.}

\item{x}{A bootstrap result to print.}
\item{parm}{Statistic names or distinct positive indices to include in an
interval table. Omission includes all statistics.}
\item{level}{Confidence level strictly between zero and one.}
\item{allow.failed}{By default, intervals are refused if any replicate
failed. Explicit \code{TRUE} permits intervals from successful replicates
and produces a warning and an attribute stating that inference is
conditional on fitting success. At least two successes are required.}
\item{...}{Printing arguments, or additional generic arguments. Unsupported
interval arguments are rejected.}
}
\details{
Repeated fits run in the current R process. Numerical kernels use the
existing OpenMP and BLAS thread controls; \code{\link{st.set.num.threads}} describes these controls.
R callbacks execute on the main R thread. Intermediate state stays in memory.

Bootstrap stream scoping requires a generator whose complete state can be
saved. Box-Muller and user-supplied generators are rejected before any RNG
state is changed. The default Inversion normal generator is supported.
The simulation covariance factorization is prepared once per run and reused
for every parametric response. Responses are generated before fitting, using the stream assigned to each
replicate. The covariance factorization is reused across responses.

The parametric method uses the fitted physical covariance matrices and
the observed fixed-effect estimate. It includes the fitted offsets and
observation precision weights. By default, random coefficients and residual
errors are regenerated independently. All model matrices, grouping indices,
factor coding, transformation calibration, known covariance matrices and
missing-response positions are retained during refitting. A subset-selection
search or an internal random-effect bootstrap test is not repeated.

When selected random coefficients are held, they are treated as fixed
plug-in constants. Unselected coefficients are generated from their marginal
fitted covariance, with both covariance axes limited to those
coefficients. This does not draw their conditional distribution given the
held coefficients. It does not propagate uncertainty in the held fitted
values. The resulting conditional target determines the interpretation of a statistic.

Residual resampling requires a univariate model with one homogeneous,
positive, independent residual variance. Heterogeneous or correlated residual
errors and multivariate responses are rejected for this method. Observed
conditional residuals are centered and scaled so that their empirical
second moment equals the fitted observation-level residual variance. They
are sampled independently with replacement and added to the full fitted
conditional mean. Constant precision weights are allowed if the resulting
residual covariance is homogeneous. This is an approximation: centering and
rescaling do not remove the shrinkage and dependence induced by fitting.
In particular, holding random effects fixed omits their between-sample
variation and can give inadequate coverage for population-level parameters.
The parametric method is the primary Gaussian mixed-model bootstrap.

With \code{series=NULL}, the statistic is an R callback evaluated on the main thread. It should use
the fitted model and explicit arguments; deliberate random draws use the
assigned replicate stream. Reproducibility requires that the callback does
not depend on mutable external state or external side effects.

The original statistic is evaluated in a separately recorded substream.
Each replicate ID has a separate \code{"L'Ecuyer-CMRG"} stream. Responses are
generated in the parent process before the repeated fits; the resulting
post-generation stream is used for a statistic that draws random numbers.
Stream assignment is determined by replicate ID. Floating-point fitting
results can differ slightly across numerical configurations.

\strong{Joint bootstrap of a series.}
Supplying \code{series} uses a parametric bootstrap which regenerates plot
responses, refits each earlier model, recomputes its adjusted means and full
covariance, and refits the second-stage model with that new covariance.
Numerical arrays, simulation roots and intermediate fitted contexts are compiled
C memory. No intermediate R fitted models or R numerical workspaces are built.
The earlier input fits and second-stage object are not changed.

The initial joint scope is univariate converged Gaussian ML/REML models,
independent first-stage errors, \code{type="parametric"}, \code{re.form=NA},
\code{statistic=stats::coef}, and no \code{statistic.args} or externally
supplied \code{responses}. Earlier provenance/order, targets, means and
within-fit covariance must match the series object. Second-stage subsets or
reordered rows, multivariate joint generation, nonzero cross-fit error covariance,
residual resampling, held random effects and custom statistics are rejected.
Explicit reference grids, named row weights, offsets, missing plot responses
and supported physical covariance constraints are retained. The target
covariance of each earlier fit must be positive definite; no ridge is added.

The matrix \eqn{L_e} defines the earlier fixed-effect targets,
\eqn{C_e=L_e\widehat{\mathrm{Var}}(\widehat\beta_e)L_e^T}, and
\eqn{\widehat m_e} includes their reference offsets. The second-stage target
signal \eqn{\theta^*} is drawn from its fitted fixed mean, random effects and
any additional estimated residual components, excluding the known inherited
covariance. For each earlier trial, the lifting matrix is
\deqn{A_e=X_e\widehat{\mathrm{Var}}(\widehat\beta_e)L_e^TC_e^{-1}}
and the plot responses are drawn with mean
\deqn{X_e\widehat\beta_e+o_e+A_e(\theta_e^*-\widehat m_e)}
and its original fitted plot-level marginal covariance. Thus target variation
from the second-stage model and first-stage plot sampling error are generated
separately; inherited sampling error is generated once. This covariance-based
lifting fixes a reproducible choice of the earlier nuisance mean when targets
use different covariate or nuisance-factor references. Earlier designs,
contrasts, transformations, masks and reference definitions remain fixed.
Every replicate re-estimates the free covariance parameters at both stages.
This is plug-in parametric generation; it is not a coverage guarantee or a
replacement for fitting a joint plot-level likelihood.

The result retains the ordinary bootstrap class and has attribute
\code{bootstrap.stage="joint series"}. The original earlier fits are supplied
explicitly when calling the function; they are not retained as backend working
objects. Each record's \code{response} is a named list of generated earlier
plot-response vectors; its \code{latent.targets} attribute records the generated
second-stage target signal. \code{record$fit$stage1} records each compiled refit's
status, covariance components, adjusted means and recomputed target covariance;
\code{record$fit$stage2} records second-stage diagnostics and covariance
components. Failed rows remain \code{NA}, with available stage diagnostics,
and are not replaced. Existing interval failure rules continue to apply.
Summary and interval tables retain the bootstrap-stage attribute.

The quantities \eqn{t_0}, \eqn{\bar t^*}, \eqn{s^*}, and \eqn{q_a} denote
the original statistic, its bootstrap mean, its sample standard deviation,
and its bootstrap quantile of order \eqn{a}, respectively. Percentile intervals use the empirical lower and
upper quantiles, with linear interpolation equivalent to R's
\code{quantile(type=7)}. Basic intervals reflect these limits around
\eqn{t_0}. Normal intervals use
\eqn{2t_0-\bar t^* \mathbin{\pm} z_{(1+level)/2}s^*}, including the
bootstrap bias correction. These are approximate intervals; nominal
coverage is not guaranteed in small samples or at variance boundaries.
BCa and studentized bootstrap intervals are not implemented by this interface.
}
\value{
An \code{st.lmm.bootstrap} object contains \code{t0}, the original named
statistic; \code{t}, the replicate-by-statistic matrix with failed rows
retained as \code{NA}; \code{records}, the full per-ID responses, streams,
statuses, warnings, errors, diagnostics and execution information; and
\code{seed}, the original model and
statistic call. Original-statistic warnings are stored separately.

\code{summary} returns the original estimate, bootstrap mean, bias, sample
standard deviation, and counts of successes and failures. \code{confint}
adds lower and upper limits, retaining statistic names and requested order.
An interval table records its type, level and whether it conditions on
successful fitting. No failed replicate is silently replaced.
}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lmm.refit}},
\code{\link{st.lmm.simulate}}, \code{\link{st.lmm.profile}}}
