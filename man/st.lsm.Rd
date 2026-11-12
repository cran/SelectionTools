\name{st.lsm}
\alias{st.lsm}
\alias{st.lsm.cld}
\alias{st.lsm.test}
\alias{st.lsm.mvcontrast}
\alias{print.st.lsm.mvcontrast}
\alias{print.st.lsm}
\alias{print.summary.st.lsm}
\alias{summary.st.lsm}
\alias{as.data.frame.st.lsm}
\alias{coef.st.lsm}
\alias{vcov.st.lsm}
\alias{confint.st.lsm}
\alias{pairs.st.lsm}
\title{Least-squares means, contrasts and fixed-effect hypotheses}
\description{Constructs reference grids and estimates marginal linear functions of an \code{st.lmm} fit. Numerical construction and inference are implemented in C. The implementation follows the reference-grid methodology of \pkg{emmeans}; that package is needed only for optional reference comparisons.}
\usage{
st.lsm(object, specs, by = NULL, at = list(), weights = "equal",
       ddf = NULL, level = .95, adjust = "none", type = "link",
       transform = NULL, trend = NULL, delta.var = NULL,
       cov.reduce = mean, cov.keep = "2", rg.limit = 10000L,
       nesting = NULL, bias.adjust = FALSE, sigma = NULL,
       regrid = NULL, max.degree = 1L, L = NULL, rhs = 0,
       joint = FALSE, side = "two.sided", test = "F", DF = NULL,
       vcov. = NULL, mvt.control = NULL, calpha = NULL,
       method = NULL, ref = NULL, wts = NULL, scores = NULL,
       reverse = FALSE, infer = NULL, inference.space = "broad")
st.lsm.cld(object, alpha = .05, adjust = "tukey",
           sort = TRUE, reversed = FALSE)
st.lsm.test(object, joint = TRUE, by = NULL, rhs = 0, ...)
st.lsm.mvcontrast(object, method = "eff", mult.name = "trait",
                 by = NULL, adjust = "sidak", rhs = 0,
                 show.ests = FALSE, ...)
\method{print}{st.lsm.mvcontrast}(x, digits = NULL, top = NULL, bottom = NULL,
  char.per.line = NULL, tests = NULL, ...)
\method{print}{st.lsm}(x, digits = NULL, top = NULL, bottom = NULL,
  char.per.line = NULL, tests = NULL, ...)
\method{print}{summary.st.lsm}(x, digits = NULL, top = NULL, bottom = NULL,
  char.per.line = NULL, tests = NULL, ...)
\method{summary}{st.lsm}(object, level = NULL, adjust = NULL,
                        type = NULL, rhs = 0, side = "two.sided",
                        mvt.control = NULL, bias.adjust = NULL,
                        sigma = NULL, DF = NULL, calpha = NULL, ...)
\method{as.data.frame}{st.lsm}(x, ...)
\method{coef}{st.lsm}(object, ..., rhs = FALSE)
\method{vcov}{st.lsm}(object, ...)
\method{confint}{st.lsm}(object, parm, level = .95, ...)
\method{pairs}{st.lsm}(x, reverse = FALSE, ...)
}
\arguments{
\item{object,x}{An \code{st.lmm} fit for \code{st.lsm}; an \code{st.lsm} object for the other functions. The multivariate contrast print method takes a \code{st.lsm.mvcontrast} result. An existing least-squares means or contrast object is also accepted by \code{st.lsm} with \code{method} to construct contrasts, with symbolic \code{L} to define hypotheses on its displayed means, or with \code{joint=TRUE} to test its linear functions.}
\item{specs}{Predictor names, or a one-sided formula such as \code{~treatment}, \code{~treatment * year}, or \code{~treatment | year}. \code{~1} requests the overall mean. This argument is absent when \code{L} or an existing \code{st.lsm} object is supplied; the forms are mutually exclusive.}
\item{by}{Conditioning variables. Contrasts retain the object's by groups when NULL; \code{character()} combines them.}
\item{at}{Named list of reference values for predictors.}
\item{inference.space}{The population described by the least-squares
means and their SEs. The default \code{"broad"} averages over the
random effects. The value \code{"narrow"} describes the specific groups
in the experiment for all fitted random terms. A one-sided random
formula or a character vector of fitted random block names selects
which groups are specific to the experiment in an intermediate space.
The value \code{~0} also gives broad inference. Intermediate and narrow
inference use reference-grid means from ML or REML fits.}
\item{weights}{\code{"equal"}, \code{"proportional"}, \code{"outer"}, \code{"cells"}, \code{"flat"}, or nonnegative numeric weights for the averaged grid cells.}
\item{ddf}{NULL retains the fit's method; otherwise \code{"Satterthwaite"}, \code{"Kenward-Roger"}, or \code{"asymptotic"} for likelihood fits. Kenward--Roger requires REML.
For a fit using \code{method="moments"}, NULL retains the fit's choice
(default \code{"Satterthwaite"}). Single-stratum functions retain exact DF;
functions combining strata use Satterthwaite DF. An explicit \code{"exact"}
requires exact inference.}
\item{level}{Confidence level, strictly between zero and one. NULL in \code{summary} retains the object's value.}
\item{adjust}{Multiplicity method: \code{"none"}, \code{"tukey"}, \code{"scheffe"}, \code{"bonferroni"}, \code{"sidak"}, \code{"holm"}, \code{"hochberg"}, \code{"BH"} (or \code{"fdr"}, \code{"FDR"}), \code{"BY"}, \code{"mvt"}, \code{"single-step"}, \code{"free"}, \code{"Shaffer"}, \code{"Westfall"}, or \code{"hommel"}. Complete pairwise contrasts default to Tukey.
The new \code{"dunnett"}, \code{"trt.vs.ctrl1"}, \code{"trt.vs.ctrlk"} and \code{"mean_chg"} families default to mvt; \code{"del.eff"} defaults to FDR. Existing other families retain their default of none; an explicit \code{adjust="mvt"} requests adjusted \code{"trt.vs.ctrl"} inference.}
\item{type}{\code{"link"} or \code{"response"}. NULL retains the stored choice.}
\item{transform}{Response transformation: NULL detects a direct \code{log}, \code{log2}, \code{log10}, \code{log1p}, or \code{sqrt} response; \code{"identity"} disables back-transformation. May be specified for a precomputed transformed response. Also accepts \code{"inverse"}, \code{"logit"}, \code{"probit"}, \code{"cauchit"}, \code{"cloglog"}, or a named list with \code{linkinv} and \code{mu.eta} functions and optional \code{name}.}
\item{trend}{Optional numeric predictor name. Estimates its local slope instead of the least-squares mean. Slopes remain on the linear predictor scale.}
\item{delta.var}{Finite-difference increment for trends; default is 0.001 times the observed predictor range. Must be nonzero.}
\item{method}{NULL constructs least-squares means from a fitted model. Otherwise constructs contrasts of an existing \code{st.lsm} object, or of the means specified by \code{specs} in the same call. Available choices are \code{"pairwise"}, \code{"revpairwise"}, \code{"trt.vs.ctrl"}, \code{"consec"}, \code{"eff"}, \code{"del.eff"}, \code{"poly"}, \code{"opoly"}, \code{"helmert"}, \code{"mean_chg"}, \code{"identity"}, \code{"tukey"}, \code{"dunnett"}, \code{"trt.vs.ctrl1"}, \code{"trt.vs.ctrlk"}, a numeric matrix with one row per contrast, or a named list of custom
coefficient vectors. Additional sample-size weighted families are
\code{"Dunnett"}, \code{"Tukey"}, \code{"Sequen"}, \code{"AVE"},
\code{"Changepoint"}, \code{"Williams"}, \code{"Marcus"},
\code{"McDermott"}, \code{"UmbrellaWilliams"} and \code{"GrandMean"}.
These capitalized names use \pkg{multcomp} contrast conventions and default
to \code{"single-step"}; existing lowercase families retain their defaults,
scaling and ordering. Their \code{wts} option supplies positive sample-size
weights, defaulting to equal weights. Each vector has one entry per mean in each by group, in displayed grid order.}
\item{ref}{Control indices or exact displayed level labels, including a vector to compare treatments with the average of several controls. NULL selects the first mean, or the last for \code{"trt.vs.ctrlk"}.}
\item{cov.reduce}{Function, named list of functions, or logical. The default takes each numeric covariate's mean. Functions may return several reference values, for example \code{range}. FALSE retains all observed values. Formula-based covariate reduction is not implemented.}
\item{cov.keep}{Names of numeric covariates to retain at their observed values, and/or a character or numeric threshold for retaining covariates with at most that many unique values (default \code{"2"}). Retained covariates are treated as factors when computing cell weights.}
\item{rg.limit}{Positive maximum number of crossed reference-grid rows, default 10000. Larger grids require a higher limit; the covariance of many requested means can require substantial memory.}
\item{nesting}{NULL detects nesting from observed memberships and hierarchical fixed terms. A named list such as \code{list(village="region")} supplies it explicitly; \code{list()} disables detection.}
\item{bias.adjust}{Controls the second-order Taylor bias correction to response summaries or regridding. NULL in summary retains the stored choice.}
\item{sigma}{Standard deviation for bias adjustment, scalar or one per summary/reference-grid row. NULL uses the square root of observation-averaged residual variance for the corresponding trait. For heterogeneous residuals or a different target distribution, an explicit sigma specifies the intended target distribution.}
\item{regrid}{NULL or \code{"none"} retains the linear-predictor grid. \code{"response"} transforms each grid prediction before averaging, propagates delta-method covariance, and places subsequent contrasts/tests on the response scale.}
\item{max.degree}{For trends, integer 1--5, default 1. Degrees above one add a \code{degree} conditioning column. For polynomial contrast families, omitting this argument selects the smaller of six and the number of means per group minus one. An explicit value selects the contrast degree.}
\item{mvt.control}{Named integration controls: \code{maxpts=131072}, \code{abseps=1e-5}, \code{releps=0}. The compiled deterministic integrator preserves R's random-number state. Summaries expose \code{mvt.error} and \code{mvt.converged}; the printer reports when the requested tolerance was not reached.}
\item{L}{Optional numeric coefficient vector or matrix, with one column per
fixed coefficient, in \code{names(coef(object))} order. Named columns are
aligned by name. Rows need not sum to zero. Estimable row-wise hypotheses require positive
standard errors. Also accepts character or
expression vectors containing linear equalities or inequalities. Symbolic
names must be coefficient names when \code{object} is an \code{st.lmm}
fit. With an existing \code{st.lsm} object, symbolic names instead identify
its displayed least-squares means or contrasts. One factor is identified by its displayed
level names; several factors are identified by displayed values joined with single
spaces in the order of the object's label columns, excluding \code{by}
columns. Nonsyntactic names require backticks. Names must be unique within
each by group. The same hypotheses are resolved separately within every
group, preserving the object's averaging, weights, covariance and offsets.
An unknown or ambiguous name is an error; no reference level is substituted.
Numeric hypotheses on an existing least-squares means object use \code{method},
not \code{L}. Exactly one of \code{method} and \code{L} is supplied, with \code{specs} absent
when supplying an existing least-squares means object.
Oppositely directed inequalities are sign-reversed into one common family
direction; the stored matrix records that orientation. Only affine arithmetic with constants, \code{+}, \code{-}, \code{*} and
\code{/} is accepted; expressions are never evaluated as R commands.}
\item{rhs}{In summaries, the null value on the linear predictor scale; omission retains the stored right-hand side. For \code{coef}, TRUE extracts the stored null vector. Otherwise a finite scalar or one null value per row of \code{L}. Symbolic
hypotheses specify their own right-hand sides, including nonzero values;
this argument is absent for them. A symbolic equality can compare one mean with a
constant, or any affine linear combination of means with another such
combination. A zero-null test requires an explicit zero right-hand side.}
\item{joint}{For \code{st.lsm}, TRUE selects a model-based joint test,
including full hypothesis-specific Satterthwaite or Kenward--Roger inference.
This also applies to symbolic hypotheses on an existing least-squares means
object, after composing them with its stored coefficient matrix. Joint tests
require equalities. Inequalities use row-wise inference and cannot be mixed
with equalities in the same family. Regridded response-scale objects retain
the existing restriction against model-based joint tests.
For \code{st.lsm.test}, the existing least-squares means joint convention is
retained for ordinary grid objects; direct coefficient objects use the
model-based test. FALSE selects row-wise inference.}
\item{test}{\code{"F"} or \code{"Chisq"} for a direct joint hypothesis.
The chi-squared test is an asymptotic Wald test with raw likelihood
covariance, without Kenward--Roger scaling. Ordinary and grouped grid joint
tests use F. A supplied \code{DF} is unavailable for a chi-squared test.}
\item{DF}{Optional positive scalar or \code{Inf} for a common Wald reference
number of degrees of freedom. Direct coefficient construction with \code{DF}
requires \code{L} and a likelihood fit and uses Wald covariance;
\code{ddf} is absent or
\code{"asymptotic"}. Summary can explicitly select a reference DF without
changing the stored object. The contrast path also accepts this summary
override. Grid construction without \code{method} requires direct \code{L}
for a constructor DF. Ordinary joint grid tests accept \code{DF} for
likelihood fits and use the Wald covariance. Moment joint tests retain their
error strata. Supplied DF are distinct from model-derived Satterthwaite DF.}
\item{vcov.}{Optional fixed-coefficient covariance matrix or a function of
the fitted object returning that matrix. Available with direct \code{L}
and likelihood fits. Dimnames are aligned, and dimensions, symmetry,
finiteness and positive semidefiniteness are checked. The matrix is not
repaired. Inference uses a supplied \code{DF}, or the normal approximation
if \code{DF} is omitted; model-derived finite-sample DF are not applied to
an arbitrary supplied covariance. The fitted object is not modified.}
\item{calpha}{Optional positive numeric critical value for confidence limits.
It overrides the interval critical value only, not the reported p-values.}
\item{mult.name}{One or more marginal-grid factors whose levels form the response vector for multivariate contrasts; default \code{"trait"}.}
\item{show.ests}{TRUE includes the individual contrast estimates. With
FALSE, the estimates table has zero rows.}
\item{alpha}{Significance threshold for letter displays.}
\item{sort}{Sort letter-display rows by mean (default TRUE). FALSE retains reference-grid row order; letter assignment still follows the direction selected by \code{reversed}.}
\item{reversed}{FALSE (default) assigns letter \code{a} starting with the smallest displayed mean; TRUE starts with the largest. Letter assignment is independent within each conditioning or trait group. When \code{sort=TRUE}, rows follow the same direction. Ties retain reference-grid order.}
\item{wts}{Optional weights for effect contrasts or positive sample-size weights for capitalized contrast families. One entry per mean in each conditioning group, in displayed order.}
\item{scores}{Optional distinct finite scores for polynomial contrasts, one per mean in each conditioning group.}
\item{infer}{For moment fits, TRUE requests intervals and tests, FALSE retains point estimates, or a logical pair selects intervals and tests separately. NULL retains the usual inference. Likelihood fits require ordinary inference.}
\item{reverse}{Reverse contrast directions; \code{pairs} reverses pairwise differences.}
\item{side}{\code{"two.sided"} (or \code{"="}), \code{">"}, or \code{"<"}. One-sided tests return matching one-sided confidence limits. Tukey and Scheffe are restricted to two-sided inference.}
\item{tests}{Printing option. NULL omits tests against zero for plain least-squares means
and retains test columns for contrasts. TRUE includes available test columns;
FALSE omits test statistics and p-values. Estimates, standard errors,
degrees of freedom and confidence intervals remain available.}
\item{digits}{Decimal places for printed estimates and standard errors (integer 1--16; default 4). Small nonzero values and p-values are shown using bounds, and t ratios exceeding 1000 in absolute value are bounded. Scientific notation is not used. Stored numbers retain full precision.}
\item{top,bottom}{Nonnegative integer numbers of leading and trailing rows to display
in each table (defaults 10 and 5). Tables with no more than their sum are shown
in full. Otherwise an explicit omitted-row count separates the two parts.
Either argument may be \code{Inf} to show every row; zero suppresses that part.}
\item{char.per.line}{Maximum output line width (default 75; integer 40--10000).
When full pairwise comparison labels exceed the available space, the printer
shortens distinguishable factor levels and gives a key for the shortened
levels. If the table still exceeds the width, the confidence limits appear
in a separate block with the same comparison labels.
The compiled printer first reduces spacing, then decimal places (to at least two,
or one if explicitly requested), then abbreviates labels with a printed key.
Probabilities and existing numerical bounds keep the requested precision.
If needed, it splits wide tables into column blocks with repeated row numbers
and the first column. Extremely wide numbers use wrapped row records. Prose and
formulas also wrap. UTF-8 widths are measured conservatively in bytes without
splitting a character. Printing does not change stored numbers or labels.}
\item{parm}{Reserved; confidence intervals are returned for every row.}
\item{...}{Optional settings forwarded by methods. Contrast options are
named arguments of \code{st.lsm}; \code{pairs} forwards them to that path.
Ordinary joint grid tests accept a denominator-DF method through \code{ddf}.
For likelihood fits with finite model-derived DF, the ordinary grid test
uses the minimum across its independent hypothesis basis. Regridded
ordinary joint tests retain the minimum stored DF within each group unless
supplied \code{DF} provide a Wald reference.
For \code{summary}, additional arguments are accepted for compatibility
with the generic and are not used.}
}
\details{

\code{st.lsm(e, method="pairwise")} compares an existing means
object, or \code{st.lsm(m, ~treatment, method="pairwise")} to construct
the means and contrasts in one call. A named list or numeric matrix in
\code{method} acts on the displayed means within each conditioning group;
\code{L} instead acts on fitted coefficients. The two forms are mutually exclusive.
\code{pairs(e)} is shorthand for the pairwise contrast path.
\code{method} must be explicit when constructing contrasts with
\code{st.lsm}; it is omitted when requesting a joint test of an existing
object. Combining \code{method} and \code{joint=TRUE} tests the newly
constructed contrasts.

Omitted contrast settings inherit the source object's confidence level,
scale, transformation, bias adjustment, sigma and DF method. The adjustment
defaults to the selected contrast family's convention, independently of
the adjustment on the means. An explicit \code{adjust="none"} disables
multiplicity adjustment. Existing grid construction defaults are unchanged.
\code{st.lsm(e, method=...)} replaces calls to the former exported
contrast constructor. All contrast options are named
arguments of \code{st.lsm}.

Printed tables use the heading Least-squares means and the column labels \code{lsmeans}, \code{DF}, \code{F}, \code{t}, \code{p} and \code{group}. Joint tests use \code{DF1} and \code{DF2}. Stored result columns use the same statistical names. Compact-letter memberships are stored in \code{.group} and displayed as \code{group}. Compact-letter displays reserve a fixed display position for each letter, even across omitted rows. Very wide letter sets use consecutive group panels. Printed alignment does not change group memberships.

Predictor-label columns that share a name with a result column receive the
suffix \code{.predictor}; a numeric suffix distinguishes further collisions.
For example, a factor named \code{p} appears as \code{p.predictor}, while
\code{p} contains test p values. Model formulas and reference-grid
specifications retain the original predictor names.


The reference grid crosses observed factor levels and explicit \code{at} values. Numeric covariates with more than two observed values are set to their arithmetic means; numeric predictors with one or two values retain those values. Means use rows with usable model predictors and at least one observed response. For a multivariate fit, each trait is reported as a separate by family. The first requested predictor varies fastest. The grid size is checked against \code{rg.limit}. Unordered factors whose labels are all finite numbers are sorted numerically. Ordered factors and explicit \code{at} sequences retain their specified order.

\subsection{Inference space and standard errors}{
The SE of a least-squares mean depends on the population it describes.
A broad mean averages over the random effects. An intermediate or narrow
mean includes predictions for specific groups in the experiment. The
reference grid includes the predictors used by the selected random terms
and uses groups present in the fitted data. The \code{at} argument can
restrict those groups. The same grid weights apply to the fixed and
random contributions, which are calculated from one fit.

The reported SE is the prediction-error SE at the fitted covariance
parameters. It includes the covariance between the fixed estimation error
and the random prediction error. For a selected random target,
\code{ddf = "Satterthwaite"} approximates DF using the sensitivity of the
target's prediction-error variance to the covariance parameters.
\code{ddf = "asymptotic"} uses infinite DF. The default for a selected
target is Satterthwaite unless the fit uses asymptotic inference.
Contrasts of the resulting means carry the same random target and its
prediction-error covariance.

Selected random effects apply to means on the linear predictor grid.
They require fitted grouping levels and prediction recipes. Direct
coefficient hypotheses, observation kernels, trends, response regridding,
Kenward--Roger DF, joint tests and Scheffe, Shaffer and Westfall
adjustments are outside the selected-target calculation.
}
Equal weights average all reference cells. Proportional weights use the joint frequencies of the averaged predictors; outer weights use the product of their marginal frequencies. Cells weights use observed cell counts within each requested mean, and flat weights average its occupied cells equally. Numeric weights follow the order of the averaged predictors in the fixed formula, with the first varying fastest. New fits save the fitted prediction recipes for \code{poly}, \code{scale}, \code{splines::ns} and \code{splines::bs}; these bases are reused at the reference grid. Formula expressions that produce factors, such as \code{factor(x)}, retain observed input levels. Older saved fits without basis metadata must be refitted or use precomputed columns. Raw matrix/list data columns still require separate named predictor columns.

Formula offsets and a scalar or data-derived \code{st.lmm(offset=)} argument
are evaluated on the reference grid and added once to the estimated mean.
Data columns referenced only by the argument offset are included in the grid
and can be set through \code{at}. Trends include the offset's derivative and
contrasts combine its grid values with their contrast coefficients. A
training-only external offset vector has no unambiguous reference-grid value;
\code{st.lsm} reports an error for that case. No arbitrary average offset is
substituted.

Nested factors include their parents in the grid. Impossible nested combinations are excluded, and averages over nested factors are computed within parents. A nested conditioning variable also retains its parents for contrasts. Detection uses the fitted data, including when \code{at} narrows the parent levels. The detected relation is available as the object's \code{nesting} attribute.

Each reported estimate is a linear function of the fitted coefficients. Estimability is checked against the fitted design; non-estimable means and their standard errors are NA. Standard errors and denominator degrees of freedom are computed for each mean or contrast. Kenward--Roger uses adjusted covariance and hypothesis-specific degrees of freedom; For likelihood fits, Satterthwaite uses covariance derivatives and observed parameter information. No asymptotic fallback is substituted when finite-sample information is unavailable. Missing Satterthwaite derivatives can be computed from a previously fitted asymptotic model without refitting its covariance parameters.

Balanced mixed designs fitted with \code{method="moments"} support the
same least-squares means interface for split-plot comparisons. For example,
\code{pairs(st.lsm(m, ~Variety), adjust="none")} compares whole-plot
treatments averaged over dates; \code{pairs(st.lsm(m, ~Date), adjust="none")}
compares dates averaged over varieties; and
\code{pairs(st.lsm(m, ~Date | Variety), adjust="none")} compares dates
within each variety. Each final contrast selects its own error variance.
Single-stratum contrasts retain the exact mean square and DF.
\code{pairs(st.lsm(m, ~Variety | Date), adjust="none")} compares varieties
within each date and combines the whole-plot and subplot errors.

By default, means and contrasts combining independent strata use
Satterthwaite DF. If their estimated variance is
\eqn{\widehat V = \sum_s a_s MS_s}, the denominator DF are
\deqn{\nu = \frac{\widehat V^2}{\sum_s (a_s MS_s)^2/\nu_s}.}
Here \eqn{MS_s} and \eqn{\nu_s} are the observed stratum mean squares
and their DF; the coefficients \eqn{a_s} are determined by the linear
function and design. Covariances use the corresponding bilinear functions.
No model refitting or likelihood estimation is needed. Nonpositive or
numerically unresolved variances remain errors; negative raw variance
components are not truncated.

\code{ddf="exact"} restricts inference to single-stratum functions.
When an exact-only least-squares mean combines strata, \code{st.lsm} retains
its point estimate and reference grid; contrasts of that grid can still
have exact tests. Exact-only inference for a contrast combining strata
raises an error. \code{infer=FALSE} requests point estimates explicitly.

For inference with balanced multivariate mixed moment fits, each mean
or contrast involves a single trait. Its variance and DF use the corresponding
diagonal entries of the stratum mean-product matrices. The off-diagonal
entries contribute to covariance between estimates for different traits,
as returned by \code{vcov}. Older saved fits without stratum geometry
retain point estimates. A new fit supplies the geometry for inference.

Joint hypotheses use \code{st.lsm.test} on the contrast object.
A common single stratum retains its exact F test. Otherwise, Satterthwaite
joint tests combine the DF of covariance-orthogonalized one-dimensional
hypotheses. Direct \code{st.lsm(..., L=L, joint=TRUE)} uses the same
calculation. Tukey adjustment requires balanced pairwise covariance and
common DF within each family; with combined strata its reference is
approximate. Bonferroni and Holm support contrast-specific DF.

Adjustments apply separately within by groups. Tukey applies to complete pairwise families; with two means it is evaluated as the exactly equivalent t test, including denominator DF below two. Scheffe uses the rank of the family's linear functions. Holm, Hochberg, BH and BY adjust p-values; their confidence intervals use Bonferroni and carry an \code{interval.adjust} attribute. Unsupported adjustment names report an error. Multivariate-t adjustment uses the actual correlation of each contrast family, including singular pairwise families. A common-factor quadrature handles Dunnett structures; other structures use conditional t probabilities and shifted, antithetic Halton integration. For agreement with \pkg{emmeans}, finite integration DF are rounded down after adding 0.25, with a minimum of one; DF above 9999 use the normal limit. Single comparisons retain their ordinary t inference. This interface remains restricted to \code{st.lmm} Gaussian models and has no plotting interface.

Ordinary response summaries inverse-transform marginal linear predictions. Bias adjustment uses the same second-order Taylor estimand as \pkg{emmeans}, not an exact lognormal correction. \code{regrid="response"} averages transformed grid predictions instead. Without regridding, standard errors use the delta method, intervals are transformed endpoint-wise, and tests remain on the linear predictor scale. Log-scale differences with zero-sum coefficients become ratios. For square-root back-transformation, negative predictions and interval endpoints map to zero, with a zero derivative for negative predictions, following \pkg{emmeans}. Reciprocal intervals crossing zero are rejected. Trends use successive finite differences of the saved design basis. Higher-degree entries divide by the degree factorial, following \code{emmeans::emtrends}; they remain on the linear predictor scale.

For likelihood fits, joint tests use an orthonormal basis for the estimable hypothesis space, its full covariance, and the smallest denominator DF across that basis, following \code{emmeans::test(joint=TRUE)}. Dependent hypotheses require a zero null. The explicit new path \code{st.lsm(e, joint=TRUE)} instead selects the full model-based joint calculation. Existing calls retain the least-squares means convention. Unavailable DF remain unavailable. Regridding retains the minimum original grid DF for each average or contrast.

Multivariate mean contrasts use the covariance rank as numerator DF. With mean marginal DF \eqn{v} and rank \eqn{r}, Hotelling's statistic is converted to \eqn{F=T^2(v-r+1)/(rv)}, with denominator DF \eqn{v-r+1}, as in \code{emmeans::mvcontrast}. Adjustments are applied across these tests within the remaining by groups.

Letter-display significance thresholds do not change the stored confidence level for mean intervals. Mean intervals use the requested adjustment (Sidak when pairwise tests use Tukey).

Letter displays use insert-and-absorb on the adjusted pairwise tests. They require all comparisons to have available inference. Letters are named in increasing mean order by default and in decreasing mean order with \code{reversed=TRUE}, so the selected extreme receives \code{a}. Overlapping groups can require several letters on a mean, such as \code{ab}; more than one mean can share \code{a}. This naming changes neither the tests nor which pairs share a letter. Shared letters indicate failure to reject a difference, not evidence of equality. Letter symbols can differ from other implementations while encoding the same comparison graph.
}
\section{Reported degrees-of-freedom method}{
For fixed-effects-only models fitted by moments, least-squares means and comparisons
use the residual degrees of freedom and report \code{exact}. The inherited
denominator-DF setting remains available in the object; the display describes
the calculation used. Estimates, standard errors, confidence limits and
p-values are unchanged by this reporting choice.
For pairwise comparisons from mixed models fitted by moments, the heading
reports \code{exact} when all comparisons have certified exact inference
from the same error stratum and share its denominator degrees of freedom.
Comparisons mixing strata retain their Satterthwaite heading. Supplied DF
values do not trigger this automatic exact label. Individual least-squares means retain
their own inference setting; an exact comparison does not imply exact
inference for each separate mean.
}
\section{Printing least-squares means}{Pairwise-comparison headings identify the variables
whose means are compared. Grouping variables remain identified in their
existing group labels. Plain least-squares means omit the t-statistics and p-values
for tests against zero from the printed display by default. Contrasts retain
their test columns. Estimates, standard errors, degrees of freedom and
confidence intervals remain displayed when available.
Tables show each row's denominator-DF method when methods differ.
Trend estimates identify their predictor, and grid results report explicit
\code{at} values. Nonzero null values and one-sided alternatives are
identified in test output. Back-transformed results identify the scale of
the tests.

The footer identifies the p-value and confidence-interval adjustments used
for the displayed results. Compact-letter displays report the comparison
adjustment and the significance threshold for the letters separately.
The confidence level of the means remains distinct from that threshold.
The print option \code{tests=TRUE} includes available test columns;
\code{tests=FALSE} omits them. The same option applies to printed summaries.
These display choices leave the fitted object, summary data frames and
\code{as.data.frame} results unchanged. For example,
\code{print(st.lsm(model, ~treatment), tests=TRUE)} includes the tests.
With \code{tests=NULL}, \code{st.options()$print.tests.for.bare.lsmeans}
controls their display for bare least-squares means. Explicit \code{tests} values take
precedence; contrasts retain tests by default.}
\section{Compact LSD and HSD summaries}{Balanced categorical analyses can
retain one or two common pairwise-difference thresholds. The named numeric
vectors \code{LSD} and \code{HSD} identify these comparison strata, for example
\code{Within env} and \code{Across env}. The compact \code{thresholds} table
contains the same one or two summaries with comparison SE, DF, family size
and alpha. No pairwise-comparison lists are retained for this feature.
When available, the table is an attribute of the summary data frame.
For a least-squares means result \code{x}, the call
\code{attr(summary(x), "thresholds")} returns this table.

For confidence level \eqn{1-\alpha}, comparison standard error \eqn{SE}
and denominator degrees of freedom \eqn{DF},
\deqn{LSD=t_{1-\alpha/2,DF}SE}
and, for a family of \eqn{k} requested means,
\deqn{HSD=q_{1-\alpha;k,DF}SE/\sqrt{2}.}
The studentized-range quantile is denoted by \eqn{q}. HSD equals LSD for
two means. The SE and DF belong to the difference, rather than either
individual mean. Families follow the conditioning factors and are separate
for different traits. A within/across label identifies a common comparison
stratum; it does not change the requested family size.

The summaries require complete, equally replicated categorical layouts,
equal replication of grouped random effects and at most two identifiable
common SE/DF/family-size combinations. Unbalanced data, more than two
combinations, unavailable quantiles or comparison strata without an
unambiguous factor label produce no stored LSD/HSD summaries and no footer.
Continuous-predictor models, point-only requests and contrast objects do
not acquire these summaries. Mixed-model thresholds retain the requested
comparison-DF approximation. The scale is the analysis scale; response
regridding uses its corresponding scale. Back-transformed presentation alone
does not transform additive thresholds.

Without displayed tests, printing uses at most two lines, for example
\code{Within env: LSD=12.3 HSD=15 (alpha=0.05)}.
There is no separate heading or comparison table. The complete footer is
omitted if its lines cannot fit the configured line width. The options
\code{st.options(list(print.LSD=FALSE))} and
\code{st.options(list(print.HSD=FALSE))} independently suppress the respective
values; both default to \code{TRUE}. The compact values remain stored when
printing is disabled. Summary confidence-level/DF overrides produce their
corresponding summary values without changing the source object. Existing
printed-test and compact-letter-display conventions are preserved.
}
\section{Direct and joint coefficient hypotheses}{
Direct matrices work with ML and REML fits and with balanced
orthogonal method-of-moments designs that provide coefficient inference.
For moment fits with explicit pooling, coefficient inference continues to use
the original unpooled fitted model and its independent error strata; pooled
hypothesis ANOVA tables are stored separately. With \code{joint=FALSE}, the result is an
\code{st.lsm} object: estimates remain \eqn{L\hat\beta}; stored \code{rhs}
and \code{side} define its default tests. Direct summaries display the tested
null alongside each estimate. Explicit summary arguments override
them. With \code{joint=TRUE}, \code{st.lsm} returns an object whose stored
\code{joint.tests} table is displayed by \code{print} and returned by
\code{summary} and \code{as.data.frame}.
For joint tests, consistent redundant rows are reduced to an independent
subset; inconsistent right-hand sides and non-estimable hypotheses are
rejected. Rank reduction preserves the original independent rows.
Kenward--Roger requires REML and uses adjusted covariance, denominator DF
and hypothesis-specific F scaling. Satterthwaite combines the moments of
orthogonalized one-dimensional tests. Existing fit information is reused;
missing covariance derivatives can be reconstructed when needed. No
alternative DF method is silently substituted.
For multivariate mixed moment fits, direct coefficient rows each concern
one trait. A joint F hypothesis uses rows from a single trait.
\code{st.lsm.test(e, by="trait")} produces separate joint tests for
the traits in a reference grid or contrast object. Fixed-only
multivariate moment fits provide trait-wise row inference and model-term
MANOVA tests.
}
\section{Simultaneous inference}{
\code{"single-step"} uses the joint normal or t distribution and one common
reference DF per family: a supplied \code{DF}, otherwise the minimum
marginal DF, without rounding. Direct-matrix \code{"mvt"} uses the same
policy. Unequal Satterthwaite DF do not define an exact multivariate t
law; this common-DF treatment is an explicitly reported approximation.
Existing grid-based \code{"mvt"} calls retain their previous DF convention.
\code{"free"} performs a correlation-based step-down adjustment.
\code{"Shaffer"} and \code{"Westfall"} use logical closure of homogeneous,
two-sided equalities, with a maximum of 20 hypotheses per family. Nonzero
nulls, offsets and one-sided logical-closure requests are rejected.
\code{"hommel"} completes the ordinary \code{p.adjust} method family.
Holm, Hochberg, Hommel, BH, BY and the step-down procedures do not imply
matching rectangular confidence intervals: the reported intervals use
Bonferroni, identified in \code{interval.adjust}. A supplied \code{calpha}
is identified separately. Integration controls and convergence information
are retained; the deterministic compiled integration does not change R's RNG.
}
\value{\code{st.lsm} returns a serializable S3 object. Objects for
least-squares means and individual contrasts contain the linear-function
matrix \code{L}, labels, covariance and inference settings. Objects for
joint tests contain the test table in \code{joint.tests}. The table is
empty for other \code{st.lsm} objects. This table is displayed by \code{print} and returned by \code{summary} and
\code{as.data.frame}. For direct joint F tests, its columns are \code{DF1},
\code{DF2}, \code{F}, \code{F.scaling}, \code{p} and \code{method}. For
direct chi-squared tests, the columns are \code{DF}, \code{Chisq}, \code{p}
and \code{method}.

\code{st.lsm.test} returns a test data frame. \code{st.lsm.mvcontrast}
returns a list with \code{estimates} and \code{tests} data frames.
Summary and confidence-interval methods for means and individual contrasts
return a \code{summary.st.lsm} data frame with estimates, SE, DF, confidence
limits, t ratios and p-values for the requested alternative. \code{coef}
and \code{vcov} return linear predictor estimates and their covariance for
means and individual contrasts. \code{st.lsm.cld} returns a
\code{summary.st.lsm} data frame with an additional \code{.group} column.}
\references{McLean, R. A., Sanders, W. L. and Stroup, W. W. (1991).
A unified approach to mixed linear models. The American Statistician, 45, 54--64.
\doi{10.1080/00031305.1991.10475767}.

Hothorn, T., Bretz, F. and Westfall, P. (2008). Simultaneous inference in general parametric models. Biometrical Journal, 50, 346--363. \doi{10.1002/bimj.200810425}.

Lenth, R. V. (2016). Least-squares means: the R package lsmeans. Journal of Statistical Software, 69(1), 1--33. \doi{10.18637/jss.v069.i01}.}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lsm-hypotheses}}}
