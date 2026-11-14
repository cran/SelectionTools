\name{st.predict}
\alias{st.predict}
\alias{print.st.predict}
\alias{print.summary.st.predict}
\alias{summary.st.predict}
\alias{as.data.frame.st.predict}
\alias{coef.st.predict}
\alias{vcov.st.predict}
\alias{confint.st.predict}
\alias{st.marker.results}
\title{Predict Random Effects and Combined Level Means}
\description{Predicts selected random effects alone or additive combinations with the fixed intercept in an \code{st.lmm} fit. A single unprefixed random-term name retains the level-mean interface. Independent random effects and supplied, possibly singular, relationship matrices are supported. Construction and numerical calculations are implemented in C.}
\usage{
st.predict(object, term = NULL, level = .95, newdata = NULL)
\method{print}{st.predict}(x, digits = NULL, top = NULL, bottom = NULL,
                          char.per.line = NULL, ...)
\method{print}{summary.st.predict}(x, digits = NULL, top = NULL,
                                  bottom = NULL, char.per.line = NULL, tests = NULL, ...)
\method{summary}{st.predict}(object, level = NULL, ...)
\method{as.data.frame}{st.predict}(x, ...)
\method{coef}{st.predict}(object, ...)
\method{vcov}{st.predict}(object, ...)
\method{confint}{st.predict}(object, parm, level = .95, ...)

st.marker.results(data.set = "default")
}
\arguments{
\item{object,x}{An \code{st.lmm} fit or the name of a marker dataset containing stored symbolic-A prediction state for \code{st.predict}; a prediction object for its methods, or a prediction summary for its print method.}
\item{newdata}{For marker prediction, the name of a target marker dataset or a supported target hybrid dataset. NULL with a marker-dataset object returns stored training conditional predictions.}
\item{term}{NULL for marker-dataset prediction, otherwise a random-term name or an additive prediction expression, for example \code{"gen"} or \code{"overlay(mother, father)"}. The selected term must be a random intercept with a constant unit basis, or a parental overlay.}
\item{level}{Pointwise coverage level, strictly between zero and one. NULL in \code{summary} retains the stored value.}
\item{digits,top,bottom,char.per.line}{The same compiled printing controls as \code{\link{st.lsm}}. Defaults show four decimal places, ten leading and five trailing rows, within 75 characters per line. Either row limit may be \code{Inf} to show all rows. Printed degrees of freedom use \code{DF}. Printing does not alter stored values.}
\item{tests}{For prediction summaries, NULL preserves available test columns,
TRUE includes them and FALSE omits them. Stored results are unchanged.}
\item{parm}{Reserved; intervals are returned for all levels.}
\item{...}{In \code{print} and \code{as.data.frame}, optional \code{level} for the summary; otherwise reserved.}

\item{data.set}{Name of a marker or supported hybrid dataset used as a source of a symbolic genomic relationship matrix in \code{st.lmm}.}
}
\details{
Additive prediction expressions select scalar random-intercept terms and
parent columns of a fitted shared overlay. For example,
\code{st.predict(fit, "0 + cross")} returns SCA effects alone, whereas
\code{st.predict(fit, "male + female + cross")} returns the intercept
plus both shared parental GCA effects and SCA. The parent names resolve
to the same overlay coefficient vector; they do not define separate effects.
Expressions use simple column names separated by \code{+}; a leading
\code{0 +} omits the intercept. Duplicate terms, unknown levels,
random slopes and ambiguous overlay parents are rejected. This interface
currently uses observed combinations in the fitted data, collapses identical
prediction targets and retains their first occurrence order. It does not
construct unobserved crosses or average additional fixed effects.
The returned \code{vcov} is the full joint prediction-error covariance
of the selected linear targets, including fixed--random and random--random
cross-covariances. Variance parameters are held at their fitted estimates.
These expressions currently require a univariate ML or REML model.


An overlay term may be prefixed with \code{0 +}, for example
\code{st.predict(fit, "0 + overlay(female, male)")}, to return parental
GCA effects alone. The unprefixed form returns intercept plus GCA.
Prediction-error covariance for the zero-intercept target uses the random
effect block only, including uncertainty due to estimation of fixed effects;
it does not add the intercept target or its cross-covariance terms.

The same interface accepts a complete overlay interaction block name, for
example \code{"0 + overlay(female, male):env"}. It returns every parental
effect within every independent grouping level, together with the full joint
prediction-error covariance in that displayed order. Group tuples and parents
use explicit length-prefixed labels, so names containing colons cannot merge
distinct targets. An \code{at()} overlay creates separate fitted blocks;
their exact names, including the stratum suffix, identify the requested blocks. An unprefixed request
adds only the fixed intercept, retaining the usual fixed-predictor reference;
it does not add an environment fixed effect automatically. Conditional
\code{predict} with the target cross rows provides complete environment-specific
cross means. Unused parents supplied through the parental relationship matrix
are retained in each grouping level.

For parental overlay predictions, printing identifies numbered parents as P1, P2, and so on
and distinguishes GCA effects alone from intercept plus shared GCA.
Stored parent identifiers, estimates, covariance and summary-column names
are unchanged.

For a series fit with only fixed environment effects, the shared genotype
term returns the equal average of fixed environment means plus its BLUP.
The full prediction-error covariance includes estimation of that average.
A separate \code{reliability} attribute (also shown in the summary) is the
variance of the genetic BLUP divided by its prior variance; it excludes
uncertainty in the average environment effect. These calculations condition
on fitted covariance parameters. This initial series prediction requires
an estimable intercept and fixed environment effects only.

When \code{object} names a marker dataset, prediction uses its stored symbolic-A fitted state and does not require a separately retained \code{st.lmm} object. \code{st.predict(fit, newdata="target")} uses the same marker-target route from the fit's snapshot. All participating sources must retain the same active fitted state for dataset-based prediction. Source changes that make this state incompatible are rejected; an existing fit retains its own snapshot. Removing a source dataset or compacting its storage does not invalidate that fit. Target datasets must still exist when requested. compiled marker datasets own their relationship matrices, effects and genotypic predictions in C memory; returned R objects are independent copies. Removing an R object does not release backend values.

Target markers are matched by identifier and encoded using the fitted source alleles, centering and relationship scale. Target frequencies do not replace the training reference. Fixed source loci contribute zero. The new marker-target route supports univariate models with an intercept-only fixed model and genomic random intercepts, including hybrid GCA and factorized SCA terms. All fitted markers are required; additional target markers are ignored. Marker order may differ. Extending the target panel therefore does not require refitting the source model: cross relationships still use the fitted marker subset and its training reference. It returns point predictions without prediction-error standard errors. Ordinary target output has \code{gen} and \code{prediction}; hybrid output retains the target performance columns and adds \code{prediction}. This mode does not modify target phenotypes. With other fixed predictors, their required target columns must be available.

The following details concern the fitted-object, named-term interface.

For an \code{overlay(mother, father)} term, the result contains one row
per parent in the common parent basis. Unless prefixed with \code{0 +},
it returns the fixed intercept plus
that parent's GCA, including their joint prediction-error covariance.
Conditional \code{predict} with parent columns in \code{newdata}
provides cross predictions.

For each level \eqn{g}, the prediction is \eqn{\widehat\mu+\widehat u_g}. Other random terms are excluded from the prediction, although they remain in the fitted model and its uncertainty calculation. For example, \code{st.predict(m, "gen")} excludes block deviations. If the fixed formula contains other predictors, only its intercept coefficient is included: numeric predictors are at zero and factor contrasts at their reference setting. This is not an average over a fixed-effect reference grid. An independently estimable fixed intercept is required.

The first implementation supports univariate Gaussian fits and random-intercept or parental overlay terms. A selected random slope, random-regression block, kernel term or multivariate response reports an error. Predictions remain on the fitted response scale, including any response transformation; there is no automatic back-transformation. All fitted levels are returned, including additional levels supplied through \code{ucov}. Factor levels retain their specified order, except that unordered numeric labels are sorted numerically, as in \code{st.lsm}. Additional nonnumeric levels follow the observed factor levels.

Genomic prediction of genotypes without phenotypes requires those genotypes
in the row and column names of the supplied relationship matrix. Their
responses may all be \code{NA}; these responses do not contribute to fitting,
but \code{st.predict} still returns the genotype predictions and their
prediction-error standard errors. For a comparison of two-step and one-step
genomic prediction requires masking every plot response of the same held-out genotypes
before either training analysis. In the two-step analysis, the field
model fits the remaining plots with fixed genotype effects and produces their
BLUEs. Only those training BLUEs supply observed responses in a genomic model
with \code{y ~ 1}, \code{random = ~gen} and \code{ucov = list(gen = A)};
the held-out means remain missing. In the one-step analysis, the genomic
model fits the same masked plot data directly, retaining the block terms.
BLUEs from a separate complete-data field model may be used as an evaluation
reference, but must not supply the training responses. Genotype identifiers determine prediction alignment; prediction row order
need not equal the data order.

Standard errors are square roots of prediction-error variances for the latent level means, including fixed-intercept uncertainty and fixed/random cross-covariance. They do not include the residual variance of a future observation. Writing \eqn{V=\operatorname{Var}(y)}, \eqn{G=\operatorname{Var}(u)}, \eqn{K=\operatorname{Cov}(y,u)}, \eqn{C=(X'V^{-1}X)^{-1}} and \eqn{L} for the repeated intercept-selection row, the prediction-error covariance is
\deqn{G-K'V^{-1}K+(L-K'V^{-1}X)C(L-K'V^{-1}X)'.}
The computation uses an estimable fixed-effect basis and does not invert \eqn{G}; zero random variance and singular relationship matrices are therefore supported.

Intervals use the normal approximation (\code{DF = Inf}) with the fitted variance parameters held fixed. They are pointwise empirical BLUP prediction-error intervals for latent means; they do not adjust for estimation of variance components. No fixed-effect t tests or p-values are attached to these random-effect predictions. The fit's Satterthwaite or Kenward--Roger setting is not used for these intervals.
}
\value{\code{st.predict} returns a serializable \code{st.predict}
object containing \code{term}, \code{levels}, \code{estimate},
\code{vcov}, \code{level} and \code{point.predictions}.
\code{coef} returns a prediction vector and \code{vcov} returns a matrix.
For a named random term, the matrix contains the full prediction-error
covariance and the point-prediction table has zero rows. The summary
is a \code{summary.st.predict} data frame with a level-label column,
\code{prediction}, \code{SE}, \code{DF}, \code{lower.CL} and
\code{upper.CL}. The level column uses the term name, or \code{level}
if that name collides with a numeric output column.

For a marker-target call, \code{summary} and \code{as.data.frame}
return the point-prediction table with the target identifiers.
A dataset call without a target gives the training conditional
predictions in the table, with one column per trait. For multiple
traits, \code{coef} stacks the columns of this table in order.
For both marker routes the prediction-error covariance is a
zero-by-zero matrix.}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lsm}}}
\examples{
data("v-lmm-rapeseed", package = "SelectionTools")
m <- st.lmm(hgt ~ env, ~geno, data = rapeseed,
            ddf = "asymptotic", random.test = "none")
p <- st.predict(m, "geno")
p
coef(p)
vcov(p)
st.predict(m, "0 + geno")
}

\section{Stored marker prediction state}{
Returns the active fitted effects and prediction information stored by a successful symbolic-A fit.

A successful top-level fit with quoted symbolic \code{"A"} or \code{"kro(a,a)"} replaces the participating source datasets' previous state. Failed or nonconverged fits do not replace it. Two parental sources of a hybrid fit share one model identity. Numeric covariance matrices and fits without symbolic covariance requests do not establish this association. Stored state is separate from the older \code{gs.return.effects} result. This accessor returns a private copy.

A versioned list containing \code{model}, \code{sources}, marker-reference \code{panels}, \code{dataset}, source hybrid \code{performance} where applicable, \code{genotype.effects}, \code{marker.effects}, and training conditional \code{predictions}. Marker effects are available for supported single-pool genomic terms; SCA is retained through relationship scores and parental kernels. The model identity is an attribute.

}
