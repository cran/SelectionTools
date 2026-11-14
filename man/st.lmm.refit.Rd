\name{st.lmm.refit}
\alias{st.lmm.refit}
\alias{st.lmm.objective}
\title{Refits a Retained Gaussian Model and Evaluates Its Likelihood}
\description{
An ML or REML model is refitted with replacement responses, observation precision
weights, or an estimation method.  The fixed-effect design, random-effect
bases, grouping indices, known covariance matrices, factor contrasts and
offset values are retained.  A separate interface evaluates the Gaussian
likelihood and its covariance-parameter gradient without fitting a model.
}
\usage{
st.lmm.refit(object, response, weights, method = NULL,
             control = list(), warm = TRUE)

st.lmm.objective(object, parameters = NULL, method = NULL,
                 gradient = TRUE)
}
\arguments{
  \item{object}{An \code{st.lmm} object fitted by ML or REML.}
  \item{response}{Optional replacement response on the original response
  scale, before subtraction of offsets. A numeric matrix has one row per
  retained observation and one column per trait. A univariate model also
  accepts a numeric vector. Values are matched by position and the
  missing-value pattern follows \code{model.response(model.frame(object))}.  Omitting this argument retains the response; explicit
  \code{NULL} is an error.}
  \item{weights}{Optional finite, strictly positive numeric precision
  weights, one per retained data row.  This argument is supported for
  univariate models.  Omitting it retains the fitted weights; explicit
  \code{NULL} removes them on the retained rows.}
  \item{method}{\code{NULL} to retain the fitted method, or \code{"ML"}
  or \code{"REML"}.}
  \item{control}{A named list overriding the retained fitting controls.
  Further details appear in \code{\link{st.lmm}}.}
  \item{warm}{Logical.  TRUE uses the fitted covariance matrices, in their
  original units, as starting values.  \code{FALSE} uses the usual
  default starts for estimated covariance blocks.  Fixed covariance
  blocks are retained.}
  \item{parameters}{A numeric vector in the coordinate system reported
  by \code{st.lmm.objective(object)}.  \code{NULL} uses the fitted
  parameters.  If names are supplied, they must match the reported names
  and order exactly.}
  \item{gradient}{Logical.  A value of \code{TRUE} includes the gradient of minus twice the
  log likelihood with respect to the reported parameters.}
}
\details{
Refitting reuses the saved model representation and does not reevaluate
the original data binding, predictor expressions, subset, or missing-value
action.  The parent fit and supplied response and weight vectors remain
unchanged.  Rows excluded when the parent model was fitted cannot be
restored by refitting.  This also applies to rows excluded by zero ordinary
regression weights.  Missing-response rows retained for prediction remain
in the refitted model; replacement responses must preserve their missing
values.  Changes to the observation layout or model specification require
\code{\link{st.lmm}} or \code{\link{update}}.

Precision weights modify only the residual covariance.  Its row and column
\eqn{i} are scaled by \eqn{1/\sqrt{w_i}}; random-effect covariances retain
their original meaning.  A response replacement is the evaluated response:
a formula such as \code{log(y) ~ x} requires new log responses.  Offsets
are subtracted internally.  The evaluated replacement response is retained
through serialization and through updates using the same response
expression and data snapshot.  An explicit new data argument or a changed
response expression starts a new response evaluation.

Replacement weights, including explicit removal by \code{weights = NULL},
are retained through serialization and updates using the same data snapshot.
They are matched to the retained row layout, including repeated rows created
by a subset.  An explicit \code{weights} argument to \code{update} overrides
these saved weights.  Replacement data after a numeric weight refit
require explicit weights for the new data.  Previously removed weights
remain removed when replacing data unless new weights are supplied.

Changing a Kenward--Roger fit to ML produces a warning and uses asymptotic
inference, because Kenward--Roger inference requires REML.  Other inference
settings are retained.  A model fitted with subset selection is refitted as
its full fitted model, with a warning; the subset search is not repeated.

Warm starts transfer covariance matrices in physical units.  Optimizer
coordinates are scaled using the current response and residual bases and
can therefore change after replacing responses or weights.  Previously
stored \code{control$parameters} are discarded.  Explicit parameters in
the new \code{control} list override the warm start and must use the new
fit's coordinate system.
For structured random-effect covariances, the direct correlation coordinate
is retained as well.  This preserves its sign even when a zero variance
makes the fitted covariance matrix insufficient to identify that sign.

The objective is minus twice the full Gaussian ML or REML log likelihood,
including its normalization constant.  Fixed effects are profiled out;
covariance scale remains an explicit parameter.  The function does not
optimize parameters, change the model, or return an external pointer or a
stateful function.  Its numeric record can be serialized or used by an
external optimizer.  Parameters outside the stated bounds are rejected.
A point satisfying the parameter bounds but failing to define a valid
marginal covariance returns \code{value = Inf}, \code{valid = FALSE} and
an \code{NA} gradient.

The coordinate table describes each variance, common variance, standard
deviation, Cholesky entry, log variance, log standard deviation, direct
correlation, or partial-correlation coordinate.  Its physical
value is the reported \code{scale} times the function in \code{inverse}
applied to the parameter: \code{identity}, \code{exp}, or \code{tanh}.
The \code{covariance} component gives the resulting covariance matrices
in the original coefficient units.  Coordinates are specific to the
retained model and must be obtained again after a response or weight refit.
Row and column zero denote a block-wide common variance or correlation.
Compound-symmetry correlation bounds depend on the number of coefficients;
AR1 correlations have bounds minus one and one.  These direct correlations
use \code{identity}, rather than a Fisher transformation.

Fixed/shared physical covariance constraints are retained by refitting and
recompiled for the response and weight scaling of the new fit. Warm and
default starts are reconciled with the retained declarations. The objective
accepts the complete raw coordinate vector as before, but a supplied vector
must satisfy every retained fixed/shared constraint. Violations are explicit
errors; arbitrary evaluation vectors are not silently projected.

The additional \code{constraints} record reports canonical physical SD and
correlation selectors, variance aliases, physical values, fixed/free/shared
status, complete and independent parameter counts, and each raw coordinate's
independent index. \code{free.parameters}, \code{free.coordinate} and
\code{free.gradient} describe the reduced optimizer coordinates. In an
objective result, both reported gradients differentiate the same minus-twice-
log-likelihood value. A shared coordinate's reduced derivative is the chain-
rule sum over its members; fixed coordinates do not contribute a free direction.
The original-scale constraints are distinct from the complete normalized
optimizer coordinates. The \code{\link{st.lmm}} help page describes the structures
and operations.
}
\value{
\code{st.lmm.refit} returns a complete \code{st.lmm} object with recomputed
estimates, covariance matrices, tests and predictions, together with a
\code{refit} record describing the replacement arguments and initialization.

\code{st.lmm.objective} returns a list containing \code{value},
\code{gradient} (or \code{NULL}), named \code{parameters}, named
\code{lower} and \code{upper} bounds, the \code{coordinates} data frame,
physical \code{covariance} matrices, \code{method}, \code{nobs}, fixed-effect
\code{rank}, a \code{valid} flag, a description of the likelihood scale,
and the physical \code{constraints} inventory and reduced-coordinate metadata.
The coordinate table identifies the covariance block, its type and
structure, coefficient row and column, transformation and physical scale.
For separable products there is one common variance coordinate and one
\code{rho.axis} coordinate for each correlated axis. Axis-specific CS bounds
are retained. Independent axes add no correlation coordinate. Refitting retains
the product layout and direct-sum strata.
For US trait-by-occasion products, the trait coordinates are
\code{trait.log.sd} and \code{trait.partial.correlation.z}, followed by
\code{rho.occasion} for AR1 or CS. Trait-coordinate row and column indices
refer to traits, not expanded trait-by-occasion coefficients. Warm starts
transfer the complete physical product matrix; cold starts preserve the
family and declared occasion grid. Missing trait cells and missing occasions
are retained under the fitted missing-data rules. Explicit covariance
constraints are not available for this US family.
}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lmm.components}},
\code{\link{update}}, \code{\link{optim}}}
