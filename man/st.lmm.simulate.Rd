\name{st.lmm.simulate}
\alias{st.lmm.simulate}
\alias{simulate.st.lmm}
\title{Simulates Gaussian linear mixed-model responses}
\description{
Generates responses from a fitted likelihood model or from a fully specified
Gaussian model without fitting it. Selected random coefficients can be held
at their fitted values while the other coefficients and residual errors are
regenerated.
}
\usage{
\method{simulate}{st.lmm}(object, nsim = 1, seed = NULL, re.form = NA,
    conditional = NULL,
    newdata = NULL, newparams = NULL, allow.new.levels = FALSE,
    offset = NULL, weights = NULL, ...)

st.lmm.simulate(fixed, random = NULL, cov = ~units, data, parameters,
    nsim = 1, seed = NULL, weights = NULL, offset = NULL, ucov = NULL,
    contrasts = NULL, traits = NULL, ...)
}
\arguments{
\item{object}{An \code{st.lmm} object fitted by ML or REML. Simulation from
the moments method is not supported.}
\item{nsim}{A positive integer giving the number of response replicates.}
\item{seed}{\code{NULL}, or a finite integer passed to \code{set.seed}.
An explicit seed is used within the call and the caller's previous random
number state, including its absence, is restored. Seeded calls reject
Box--Muller and user-supplied generators because their complete state is not
stored in \code{.Random.seed}. With \code{NULL}, the call advances the
caller's current stream without restoring it; these generators are then
allowed, but replay from the returned seed alone is not guaranteed.}
\item{re.form}{Random coefficients to hold at their fitted values. The default
\code{NA}, \code{~0}, \code{~-1}, or \code{character(0)} regenerates all
random coefficients. \code{NULL} holds all random coefficients. A character
vector selects exact complete covariance-block names. A one-sided random
formula can select fitted coefficients within a block: for example,
\code{~1|group} holds only the random intercept from a fitted
\code{~1+x|group} term. Formulas combining several selected terms are
supported. Names and coefficient definitions must match the fitted model.}
\item{conditional}{Optional string convenience policy. \code{"none"} is
exactly equivalent to \code{re.form=NA} and regenerates all random
coefficients; \code{"all"} is exactly equivalent to \code{re.form=NULL} and
holds all fitted random coefficients. \code{NULL} preserves the established
\code{re.form} argument and default.}
\item{newdata}{A data frame or named list containing target model variables,
or \code{NULL} for the retained training rows. Its response column is not
used. Fixed predictors, needed random predictors, residual strata and any
residual occasion and grouping variables must be available and nonmissing.
Numeric model variables must be finite. Learned transformation and contrast
definitions are retained from the fit.}
\item{newparams}{Optional physical parameters, in the layout described below.
Omitted top-level entries use fitted values. No covariance optimization is
performed.}
\item{allow.new.levels}{Whether new levels of an independent random grouping
factor are allowed. They have zero fitted random means and their specified
random covariance when regenerated. A named known covariance matrix must
already contain the target levels, or the fitted hybrid covariance must
define their parental relationships.}
\item{fixed}{For \code{st.lmm.simulate}, a one-sided fixed-effects formula.
There is no observed response and no response expression is evaluated.
\code{~0} specifies a model with no fixed coefficients.}
\item{random, cov}{Random-effect and residual covariance formulas, using the
same syntax and covariance families as \code{st.lmm}. Homogeneous factor
products and \code{dsum(~covariance | stratum, group=Subject)} use the
fitted axis order and independent-block structure. Supplied covariance
matrices must obey the product shape; arbitrary positive-semidefinite
matrices are not interchangeable with that family. Univariate parental
overlays support ordinary grouping interactions and selected \code{at()}
strata, retaining shared parent identities, self-cross loading two, and a
common parental relationship matrix. Independent grouping levels use separate
parental vectors; an \code{at()} block retains its own variance. Overlay
coefficient formulas, correlated coefficient products and multivariate
overlays remain unsupported.}
\item{data}{A named data frame or list of model variables. The same bitpacked
SelectionTools marker and hybrid data-set inputs accepted by \code{st.lmm}
can also supply these variables and known covariances. supported names
must be quoted. Supported specifications include \code{ucov=list(gen="A")} and the hybrid product
\code{"kro(a,a)"} for symbolic covariances; bare matrix names evaluate R
objects. The \code{random} interface also accepts \code{||}, equivalent
to a diagonal coefficient covariance with \code{|}. No observed response
is required.}
\item{parameters}{Complete physical parameters for \code{st.lmm.simulate}.
Both the named fixed-coefficient vector and all random and residual
covariance blocks must be supplied. Further details appear in Details.}
\item{weights}{For \code{st.lmm.simulate}, one finite, strictly positive
residual precision weight per row, evaluated in \code{data} and then the
formula environment. For fitted simulation, an explicit override of
future residual precision requires \code{newdata}; \code{weights=NULL}
removes target weighting. Training-row simulation retains the fitted
weights. A data-derived fitted weight expression can be reused for new
data; an external vector or a retained-weight refit requires explicit
target weights. Weights currently require a univariate response.}
\item{offset}{For \code{st.lmm.simulate}, a numeric scalar or one value per
row, evaluated in \code{data} and then the formula environment. For fitted
simulation, an explicitly supplied value replaces only the argument-offset
part; \code{NULL} removes that part. Formula offsets are added once.
An omitted value reuses the fitted offset or its stored new-data recipe.
An external training offset vector requires an explicit override for new
data. Offsets currently require a univariate response.}
\item{ucov}{Named known covariance matrices or supported SelectionTools
covariance requests, as for \code{st.lmm}.}
\item{contrasts}{Optional factor contrasts, as for \code{st.lmm}.}
\item{traits}{\code{NULL} for a univariate response, or unique, nonmissing,
nonempty trait names. Multiple names define a multivariate response without
requiring an observed response matrix.}
\item{...}{Additional arguments are rejected.
In particular, \code{st.lmm.simulate} has no fitted random coefficients to
hold and does not accept \code{re.form}.}
}
\details{
Multivariate residual products such as
\code{cov=~us(trait):ar1(Visit) | Subject} use one unstructured trait
covariance and one unit-diagonal occasion correlation matrix. Homogeneous
CS and homdiag occasions are also supported. The required input is the full physical
trait-by-occasion covariance matrix, with occasion varying fastest inside
trait and matching coefficient names. The matrix must obey the requested
product form; a general positive-semidefinite matrix is insufficient.
Declared occasion levels, including absent levels, retain their lag spacing.
Singular positive-semidefinite products are supported for specified
simulation. Exact linear relations in such draws can differ by ordinary
floating-point eigendecomposition roundoff. Further details appear in \code{\link{st.lmm}} for
the residual-only syntax and current composition restrictions.

The simulation model is \eqn{y=X\beta+o+\sum_k Z_k u_k+e}. Fixed coefficients
and covariance parameters are treated as specified values. No uncertainty
in their estimation is drawn. Unless selected by \code{re.form}, random
coefficients have zero mean and their specified block covariance; residual
errors have zero mean and the specified residual covariance, including
precision weights, residual heterogeneity and the structured residual
families described in \code{\link{st.lmm}}. With repeated-observation
residuals, occasion levels and their order are retained from the fitted
model, even when \code{newdata} contains only a subset of them. New levels
of the residual grouping factor are allowed; new occasion levels and residual
strata without fitted covariance parameters are errors. AR(1) gaps retain
their fitted factor-level distances.

The meaning of \code{re.form} differs from prediction: it selects random
coefficients to \emph{hold}, rather than components to include in a predicted
mean. Held coefficients are plug-in constants equal to their fitted values.
The remaining coefficients are regenerated from their \emph{marginal}
covariance. If an intercept is held while a correlated slope is regenerated,
the slope retains its marginal variance, with no intercept--slope covariance
in the simulated response. This is not a conditional Gaussian draw given
the held coefficients and is not a draw from their posterior distribution.

Physical parameters have the following layout:
\preformatted{
list(fixed = c("(Intercept)" = 2, x = 0.5),
     covariance = list(random = list(group = S),
                       residual = list(units = E)))
}
The names of \code{fixed} identify every fixed model-matrix column exactly
once. In multivariate models these names have the form
\code{trait:coefficient}, in trait-major order. Supplied values must be
finite, including positions aliased in a fitted model; \code{numeric(0)}
is used when the model has no fixed coefficients. Every supplied covariance
section names all blocks in that section. Block names match the fitted
model, and each block is a matrix whose row and column names match its
coefficient axes. A scalar variance is therefore a named one-by-one matrix,
not a standard deviation. Names can be reordered and are aligned explicitly.
The \code{fixed} and \code{covariance} entries of
\code{st.lmm.components(object)} provide the fitted layout. For simple
univariate intercept blocks, the matrix axis name is \code{"(Intercept)"}.
For multivariate intercept blocks, the axes have the form
\code{trait:(Intercept)}. Known
covariance between group levels is supplied through \code{ucov} separately
from these coefficient covariance matrices.

Covariance matrices must be finite, symmetric, positive semidefinite, and
consistent with the declared covariance family. \code{st.lmm.simulate}
validates a fully specified model without evaluating a likelihood or solving
for fitted coefficients. Zero variances and singular positive semidefinite
covariances are therefore allowed, including deterministic responses when
the complete covariance is zero. Fitted simulation also bypasses the observed
data solve when \code{newparams$fixed} supplies a known mean and no random
coefficients are held. Otherwise, supplied covariance parameters require a
positive definite observed covariance: omitted fixed coefficients are
estimated by GLS, and held random coefficients are recomputed as BLUPs at
the supplied parameters. Supplying known fixed coefficients recomputes these
BLUPs relative to that mean. No covariance optimizer is used.

All retained response cells are simulated, including cells whose original
responses were missing. Removed rows from an explicit \code{na.omit} or
\code{na.exclude} fit are not reinserted; \code{newdata} specifies a
different target layout. The missing-response mask is returned so a
subsequent refit or bootstrap can reinstate the original observation pattern.
Simulation needs the full covariance of the requested response cells, so
its memory use grows quadratically with their number even when the fitted
model used a sparse solver. New residual strata without fitted covariances,
and new targets for direct observation kernels without a target-kernel
definition, are rejected.
}
\value{
A list of \code{nsim} matrices named \code{sim_1}, \code{sim_2}, and so on.
Each matrix has one row per observation and one named column per trait,
including models with a single trait. For a simulation from a specification
with no response name, the column is named \code{response}. Row names follow
the target data. Attributes include \code{seed}, the actual initial R random
number state (this records \code{.Random.seed}, not any generator state
external to it); \code{missing.response}, the original response-cell mask in the same matrix layout;
\code{random.held}, a named indication of blocks with any held coefficients;
\code{random.held.coefficients}, the exact named coefficient masks;
\code{fixed.known}, the fixed-parameter interpretation; and
\code{conditioning}, a description of the simulation convention.
\code{input.rows} identifies retained original training rows when available.
}
\seealso{
\code{\link{st.lmm}}, \code{\link{st.lmm.components}},
\code{\link{st.lmm.refit}}, \code{\link[stats]{simulate}}.
}
\keyword{models}
