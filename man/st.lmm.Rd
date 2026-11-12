\name{st.lmm}
\alias{st.lmm}
\alias{st.options}
\alias{recover_data.st.lmm}
\alias{emm_basis.st.lmm}
\alias{.st_lmm_emm_DF}
\alias{modelparm.st.lmm}
\alias{print.st.lmm}
\alias{print.st.lmm.result}
\alias{Ops.st.lmm.result}
\alias{Math.st.lmm.result}
\alias{Summary.st.lmm.result}
\alias{as.data.frame.st.lmm.result}
\alias{as.matrix.st.lmm.result}
\alias{as.vector.st.lmm.result}
\alias{as.list.st.lmm.result}
\alias{[.st.lmm.result}
\alias{[[.st.lmm.result}
\alias{t.st.lmm.result}
\alias{AIC.st.lmm}
\alias{BIC.st.lmm}
\alias{confint.st.lmm}
\alias{summary.st.lmm}
\alias{print.st.lmm.selection.summary}
\alias{plot.st.lmm}
\alias{plot.st.lmm.selection.summary}
\alias{anova.st.lmm}
\alias{coef.st.lmm}
\alias{vcov.st.lmm}
\alias{fitted.st.lmm}
\alias{residuals.st.lmm}
\alias{predict.st.lmm}
\alias{logLik.st.lmm}
\alias{nobs.st.lmm}
\alias{df.residual.st.lmm}
\alias{formula.st.lmm}
\alias{model.matrix.st.lmm}
\title{Gaussian Linear Mixed Models}
\description{Fits univariate and multivariate Gaussian mixed models with formula-based specification, covariance estimation and inference.}
\section{Model data display}{
The model description begins with \code{Data:}, followed by \code{Fixed:}.
The data argument recorded in the fitted call identifies the input without
evaluating it again. A fit without a data argument reports the formula
environment, or a series of earlier fits when those provide the model data.
An embedded data snapshot reports stored data. Earlier objects
without a recorded call report that the source was not recorded.
}
\section{Fixed-effect table display}{
The response-label column \code{trait} is retained in stored fixed-effect
test tables. Model printing, \code{anova} printing and extracted test-table
printing omit this redundant column when the table contains only one trait.
Tables containing several traits retain it. This display choice does not
change the stored results or their numerical values.
}
\usage{
st.options(options = NULL)
st.lmm(fixed, random = NULL, cov = ~units, data,
       method = "REML", ddf = "Satterthwaite", type = 3L,
       random.test = "auto", control = list(), ucov = NULL,
       prediction = "NAonly", selection = NULL, weights = NULL,
       contrasts = NULL, subset = NULL, na.action = NULL, offset = NULL,
       rcov = NULL, series = NULL, constraints = NULL)
\method{print}{st.lmm}(x, digits = NULL, top = NULL, bottom = NULL,
  char.per.line = NULL, ...)
\method{print}{st.lmm.result}(x, digits = NULL, top = NULL, bottom = NULL,
  char.per.line = NULL, ...)
\method{summary}{st.lmm}(object, selection = TRUE, ...)
\method{print}{st.lmm.selection.summary}(x, digits = NULL, top = NULL,
  bottom = NULL, char.per.line = NULL, ...)
\method{plot}{st.lmm}(x, scale = NULL, labels = NULL, main = NULL,
  col = NULL, which = NULL, type = "conditional", id.n = 3L, ...)
\method{plot}{st.lmm.selection.summary}(x, scale = NULL, labels = NULL,
  main = NULL, col = NULL, ...)
\method{anova}{st.lmm}(object, effects = "fixed", ...,
                       refit = TRUE, test = "Chisq", model.names = NULL)
\method{coef}{st.lmm}(object, id = NULL, vcov = FALSE, ...)
\method{vcov}{st.lmm}(object, adjusted = NULL, id = NULL, ...)
\method{fitted}{st.lmm}(object, ..., type = "response")
\method{residuals}{st.lmm}(object, type = "response", ...)
\method{predict}{st.lmm}(object, newdata = NULL, type = "marginal",
                        se.fit = FALSE, allow.new.levels = FALSE,
                        new.levels = NULL,
                        id = NULL, offset = NULL, re.form = NULL,
                        newparams = NULL, ...)
\method{logLik}{st.lmm}(object, ...)
\method{AIC}{st.lmm}(object, ..., k = 2)
\method{BIC}{st.lmm}(object, ...)
\method{confint}{st.lmm}(object, parm, level = .95, method = "Wald", ...)
\method{nobs}{st.lmm}(object, ...)
\method{df.residual}{st.lmm}(object, ...)
\method{formula}{st.lmm}(x, id = NULL, ...)
\method{model.matrix}{st.lmm}(object, ...)
}
\arguments{
\item{options}{A named list of shared printing options, or \code{NULL} for the current values. All supplied entries are validated before any value changes.}
\item{constraints}{NULL, or a named list with \code{fixed} and/or
\code{shared} components imposing physical covariance-parameter constraints.
\code{fixed} is a named numeric vector; \code{shared} is a list of character
vectors of parameter selectors. Further details appear in \sQuote{Fixed and shared covariance parameters}
for the selector inventory, units, supported structures and inference limits.}
\item{series}{NULL, or \code{list("env", fit1, fit2, ...)} containing at
least two independent earlier ML/REML fits. The first string
names the series factor; the remaining elements are ordinary fitted
\code{st.lmm} objects with the same simple response name, or the same
set of distinct response names in \code{cbind(y1,y2,...)}. Trait order
may differ between earlier fits. Exactly one
fixed factor shared with the next-stage model identifies the adjusted
means. Additional predictors in the next-stage fixed/random formulas
are retained when they are nonmissing and constant within every trial. Further details appear in \sQuote{Series of earlier fits}.}
\item{rcov}{NULL, or a known residual covariance matrix on the absolute
response scale. Supported for univariate and multivariate ML and REML fits
without precision weights. \code{cov=NULL} fits exactly this residual covariance;
the default \code{cov=~units} adds an estimated independent residual
variance. Further details appear in \sQuote{Known residual covariance}.}
\item{fixed}{Two-sided fixed-effects formula. A matrix response such as \code{cbind(y1,y2)} creates a separate fixed-effect coefficient vector for each trait.}
\item{random}{One-sided grouping formula, a random-regression formula, or a list of these. For example, \code{~gen + env + gen:env}, \code{~Time | Chick}, \code{~diag(1+Time) | Chick}, or \code{~1+Time || Chick}. The double-bar form specifies independent expanded random coefficients and also works in additive bar terms and formula lists. Univariate random-coefficient models also support \code{homdiag()}, \code{cs()} and \code{ar1()}; further details appear in \sQuote{Random-coefficient covariance structures}.}
\item{cov}{One-sided residual covariance formula, or a list of such formulas;
default \code{~units}. Terms containing \code{units} describe independent
data rows, with possible residual covariance between traits on the same row.
For a multivariate response, \code{~units} estimates an unstructured trait
covariance matrix. \code{~diag(trait):units} specifies separate trait
variances with zero residual covariances, \code{~at(loc):units} specifies
an unstructured trait covariance matrix per location, and
\code{~diag(trait):at(loc):units} specifies a separate variance per trait
and location. Structured trait residuals also support \code{homdiag(trait)},
\code{cs(trait, hom=TRUE)} and \code{ar1(trait)}, combined with \code{units}.
For multivariate repeated records, \code{~us(trait):ar1(Visit) | Subject}
combines one unstructured trait matrix with AR1 occasion correlation; the
occasion factor also supports homogeneous CS and homdiag. Further details appear in
\sQuote{Multivariate repeated records}. For repeated observations of a
univariate response, a full factor indicator basis and an independent
group are required, for example
\code{~ar1(0 + Visit) | Subject}. This form replaces \code{units} and
allows residual covariance between different data rows in the same group.
The declared levels of the factor \code{Visit} determine the occasion order.
The supported structures and restrictions are described in
\sQuote{Covariance formula syntax}.}
\item{data}{Named data frame or list containing the model variables, or the name of a SelectionTools marker data set. A stored marker or hybrid data set must be supplied by its quoted name, such as \code{data = "maize"}. Unquoted \code{data = maize} evaluates the R object \code{maize}; an undefined name is an error. Marker datasets supply \code{y} (stored performance values) and \code{gen} (individual identifiers).}
\item{method}{\code{"REML"} (default), \code{"ML"}, or \code{"moments"}.
The method of moments uses least squares for fixed-effects models and
unconstrained moment variance and covariance estimates for supported balanced
mixed designs, including multivariate responses. With subset comparison,
\code{"REML"} still fits the specified full model by REML, but candidate
fixed-effect subsets are compared using ML. This distinction is reported in
the output; further details appear in \sQuote{Fixed-effect subset comparison}.
For \code{confint}, \code{method="Wald"} retains the default fixed-coefficient
interval; \code{method="profile"} constructs ML likelihood-profile intervals
for selected physical parameters; further details appear in \code{\link{st.lmm.profile}}.}
\item{ddf}{\code{"Satterthwaite"} (default), \code{"Kenward-Roger"}, or \code{"asymptotic"} for likelihood fits. Kenward--Roger requires REML and a supported linear covariance structure; further details appear in \sQuote{Random-coefficient covariance structures}.
With \code{"Satterthwaite"}, fixed-factor ANOVA terms in complete balanced orthogonal categorical designs use an exact expected-mean-square ratio where available and are labelled \code{"exact F"}. Other terms retain Satterthwaite inference. Explicit \code{"Kenward-Roger"} and \code{"asymptotic"} choices retain their requested reference. This automatic ANOVA choice does not change least-squares means or arbitrary coefficient-hypothesis inference.
For \code{method="moments"}, \code{"Satterthwaite"} is the default:
single-stratum linear functions retain exact DF, and functions combining
strata use Satterthwaite DF. \code{"exact"} requires exact inference.
The model's design-specific ANOVA F tests retain their exact denominators.}
\item{type}{Fixed-effect test type: 1, 2 or 3 (default). In \code{predict}, \code{"marginal"} uses fixed effects and offsets, \code{"conditional"} adds selected fitted random effects, and \code{"random"} returns only selected random effects.
Diagnostic plots support \code{"conditional"} and \code{"marginal"}.
For \code{fitted} and \code{residuals}, \code{"response"} uses the fitted
value or residual for the method's default target. This target is the
conditional mean for likelihood mixed models, the ANOVA projection for
mixed moments fits, and the marginal mean for models without random
effects. \code{"marginal"} selects the fixed mean, \code{"conditional"}
selects the mean including predicted random effects when available, and
\code{"projection"} selects the ANOVA projection for moments fits.
Likelihood fits also provide \code{residuals(type = "pearson")}. Further
details appear in \code{\link{st.lmm.diagnostics}}.}
\item{random.test}{\code{"auto"} (default) uses exact F-tests when a complete balanced orthogonal categorical design provides one denominator mean square for every random term, with an intercept, independent random intercepts and homogeneous independent errors. These automatic exact tests apply to unweighted models without offsets or supplied covariance matrices. Otherwise it uses \code{"LRT"}. The explicit choice \code{"F"} requires those exact tests; an unavailable exact denominator is an error. These tests use the same balanced projection and expected-mean-square calculations as the moments path, without changing the ML/REML fit or fitting a moments model. The exact-test table reports \code{DF1}, \code{DF2}, \code{F} and \code{p}, with \code{method="exact F"}. \code{"LRT"} explicitly tests each whole random term using an approximate chi-squared likelihood-ratio test, with one reduced-model fit per term. \code{"bootstrap"} requests a parametric bootstrap under the model without that term; \code{"none"} skips random-effect tests and their printed section.
For \code{method="moments"}, the default is \code{"F"} for univariate
responses and \code{"none"} for multivariate responses. The supported
choices are \code{"F"} and \code{"none"}. Multivariate \code{"F"} requests
separate tests by trait. An available single denominator mean square gives
an exact F-test. The general path can also use quasi-F tests based on
mean-square combinations when \code{control$quasi.f} is enabled.}
\item{control}{Named list of fitting controls: \code{optimizer="auto"},
\code{ai.maxit} (default 50), \code{profile.scale=TRUE}, \code{maxit}
(default 500), \code{factr}, \code{pgtol}, \code{lmm}, \code{trace},
\code{REPORT}; \code{bootstrap} sets the replicate count (default 999).
For a multivariate mixed model fitted by the method of moments,
\code{print.mean.products = TRUE} includes the observed mean-product
matrices (SSCP divided by DF) in the printed output. The default is
\code{FALSE}.
For an earlier-fit series, \code{series.id="plot.id"} names a common
observation-ID column and checks that observations contributing to the
first-stage fits have unique, disjoint IDs. These must identify observations
globally, not merely enumerate rows separately in each trial. Further details appear in
\sQuote{Series of earlier fits}.
The compiled convenience policies are \code{scale="none"},
\code{singular="current"}, \code{singular.tol=1e-4},
\code{boundary="current"}, and \code{start="auto"}. Scale choices are
\code{"none"}, \code{"check"}, and \code{"auto"}; all retain the engine's
internal numerical scaling, while \code{"check"} warns for a scale ratio
above 1000. Singularity choices are \code{"current"}, \code{"ignore"},
\code{"warning"}, and \code{"error"}. Boundary choices are
\code{"current"}, \code{"accept"}, \code{"restart"}, and \code{"error"}.
Start choices are \code{"auto"}, \code{"default"},
\code{"zero.correlation"}, and \code{"moments"}.
The automatic and moment policies use a projected moment initializer for
eligible univariate GCA-only models with two crossed scalar relationship terms,
one scalar independent residual, constant coefficient loadings, at least 12000
observed responses, and
\code{backend="auto"} or \code{backend="sparse"}. Automatic selection also
requires the combined relationship border to be smaller than the observed
sample. Other combinations retain the existing feasible initializer, including
an explicit \code{start="moments"} request. Supplied covariance starts and
parameter evaluations retain their values. Diagnostics
\code{start.moments.status} and \code{start.moments.rows} describe the path used.
This initialization policy is separate from \code{method="moments"}.
\code{joint=TRUE} additionally stores joint fixed/random prediction-error
covariance, including with \code{ddf="asymptotic"}.
\code{backend="auto"} retains complete spectral product shortcuts. For
eligible incomplete scalar GCA/SCA products with 6000--12000 observed responses
and at least 16 configured threads, it prefers the dense calculation. Other
models retain the spectral, sparse mixed-model, covariance-block and dense
eligibility rules. The thread count supplied by \code{st.set.num.threads}
is preserved. \code{backend="dense"} retains that reference calculation.
\code{backend="spectral"} requires an eligible spectral product model.
\code{eigen.solver="auto"} prefers DSYEVD within its workspace limit;
\code{"dsyevd"} and \code{"dsyevr"} select the driver explicitly. \code{backend="sparse"} or \code{backend="block"} forces the respective route without its size heuristic; unsupported structures are rejected. A singular base covariance falls back to the dense calculation.}
\item{ucov}{NULL, or a named list of known covariance matrices for random grouping terms, for example \code{list(gen = A)} for a genomic relationship matrix. With a stored marker dataset, the quoted string \code{"A"} requests its internally constructed VanRaden matrix; a stored hybrid interaction uses \code{"kro(a,a)"}. Unquoted \code{A} evaluates the actual R matrix. Prepared lists are supported and described in the Details section. Each entry must be a finite, symmetric positive-semidefinite square numeric matrix with unique nonempty row and column names identifying the same levels. Matrices are matched by names, including when column order differs from row order. Every non-omitted grouping level must be covered. Extra matrix levels are retained for prediction. Unmatched or duplicate list names are errors. Without an entry, the corresponding grouping levels remain independent.}
\item{x,object}{An \code{st.lmm} object, an
\code{st.lmm.selection.summary} object for its print and plot methods, or an
extracted result for the \code{st.lmm.result} print method.}
\item{parm}{Coefficient names or indices for confidence intervals; omitted for all coefficients.}
\item{level}{Confidence level for fixed-coefficient intervals.}
\item{model.names}{Optional unique nonempty names for compared models. NULL uses the supplied model expressions where concise and distinct, with numbered labels as a fallback.}
\item{refit}{For multi-model \code{anova}, TRUE refits REML models by ML. Further details appear in \code{\link{st.lmm.compare}}.}
\item{test}{For \code{effects = "manova"}, an omitted argument or NULL
returns all stored methods. The available method names are \code{"Pillai"},
\code{"Wilks"}, \code{"Hotelling-Lawley"} and \code{"Roy"}.
For multi-model \code{anova}, \code{"Chisq"} (default)
requests likelihood ratio inference where nesting is established.
For fixed-effect models fitted by the method of moments, the default
and \code{"F"} request exact F tests. The choice \code{"none"} omits
tests in either comparison.}
\item{k}{Penalty per parameter in AIC; default 2.}
\item{prediction}{For a SelectionTools hybrid data set, \code{"NAonly"} (default) predicts missing responses in the retained performance frame and \code{"Full"} predicts every combination of the two parental pools. The stored table and parental identifiers are available through \code{\link{st.lmm.components}(object, "predictions")}. Further details appear in \code{\link{st.hybrid.data.set}} for prediction coverage and covariate requirements. This option does not change non-hybrid fits.}
\item{selection}{For \code{st.lmm}, NULL leaves subset comparison disabled;
a named list enables fixed-effect subset comparison. The list controls
are described in \sQuote{Fixed-effect subset comparison}.
For \code{summary}, a logical value: TRUE (default) returns the subset
comparison summary when the fitted object contains a search, while FALSE
returns the ordinary summary of the specified full model.}
\item{weights}{NULL or a numeric vector of observation precisions,
evaluated in \code{data} and then the calling environment. With homogeneous
independent residuals, positive weight \eqn{w_i} gives residual variance
\eqn{\sigma^2/w_i}. With random effects or structured residuals, weights
require a univariate response and \code{"ML"} or \code{"REML"}; every
retained row, including a missing-response prediction row, must have a finite,
strictly positive weight. Fixed-effects-only regression with
\code{cov=~units} also supports \code{"moments"} and preserves its existing
exclusion of zero and non-finite weights. Finite negative weights are errors.
These are precision weights, not observation replication counts.}
\item{contrasts}{NULL or a named list assigning contrast matrices or
contrast functions to factor expressions. An entry may be a numeric matrix,
a function, or a function name such as \code{"contr.treatment"},
\code{"contr.sum"}, \code{"contr.helmert"}, \code{"contr.poly"} or
\code{"contr.SAS"}. It overrides the corresponding factor's contrasts
attribute. Otherwise that attribute is used, followed by
\code{getOption("contrasts")}; ordered factors use its second entry.
Matrix rows correspond to the input factor levels in order. Matrix values
must be finite; supplied column names must be nonempty. The fitted coding
is retained for prediction and subsequent updates.}
\item{subset}{Optional row selection, evaluated in the input data and calling
environment. A logical scalar or a logical vector with one value per input
row is accepted. Integer indices may be positive, including repetitions, or
negative; zeros are ignored and positive and negative indices cannot be mixed.
Missing, fractional and out-of-range indices are errors.}
\item{na.action}{NULL preserves the package's missing-response prediction
behavior, independently of the global \code{na.action} option. Alternatives
are the standard functions \code{na.omit}, \code{na.exclude}, \code{na.fail}
and \code{na.pass}, or their names as strings. Custom functions are
not accepted. The explicit policies and row mapping are described below.}
\item{offset}{For fitting, an optional known contribution to a univariate
response mean, evaluated in the input data and calling environment. It is a
numeric vector with one value per input row, or a scalar repeated over rows.
It is added to any fixed-formula offsets. For prediction, an explicitly
supplied value replaces only this argument-offset component; an explicit
NULL removes it. Prediction offsets have one value per new-data row, or per
retained model-frame row when \code{newdata=NULL}, or may be scalar.
Omitting the prediction argument uses the fitted offset recipe when possible.}
\item{id}{NULL (default) refers to the specified full model. A subset model
identifier selects a retained candidate for coefficients, covariance, formulas,
predictions or Wald intervals. Identifiers correspond to rows of the subset
summary. A single identifier or a vector of identifiers returns a list named
by identifier, in the requested order. For \code{confint}, \code{id} is
supplied through \code{...} with \code{method = "Wald"}.}
\item{vcov}{For candidate extraction with \code{coef} and \code{id}, attach
the coefficient covariance matrix when TRUE. The default is FALSE.}
\item{scale}{Criterion used to order the subset inclusion plot: \code{"rss"},
\code{"r2"}, \code{"adjr2"}, \code{"Cp"}, \code{"BICm"},
\code{"AIC"}, \code{"BIC"} or \code{"logLik"}. NULL uses \code{"BICm"}
when available and \code{"BIC"} otherwise.}
\item{labels}{Optional labels for the columns of the subset inclusion plot,
one per inclusion column, including the intercept when fitted.
For a diagnostic plot, one character label per retained response cell,
in response-matrix column order. Missing labels are not annotated.}
\item{main}{Optional plot title. For several diagnostic panels, either one
title or one title per panel.}
\item{col}{Optional colour vector for the subset inclusion plot.
For diagnostic plots, one colour or one per retained response cell.}
\item{which}{NULL retains the subset inclusion plot when the model
contains subset results and otherwise requests four diagnostic panels.
\code{"selection"} explicitly requests the subset plot. Diagnostic
selectors are \code{"residuals"}, \code{"qq"}, \code{"scale-location"}
and \code{"leverage"}, or their indices 1 through 4. Further details appear in
\code{\link{st.lmm.diagnostics}} for their definitions.}
\item{id.n}{Maximum number of diagnostic observations to annotate, chosen
by absolute standardized residual. Zero suppresses annotations.
This argument does not change subset inclusion plots.}
\item{adjusted}{For \code{vcov}, NULL uses adjusted covariance when the fit requested Kenward--Roger. TRUE requests that adjustment explicitly; FALSE returns the model-based covariance.}
\item{digits}{Number of decimal places for printed statistics and p-values (integer 1--16). \code{NULL} uses \code{st.options()$print.digits}, initially 4. Degrees of freedom use at most two decimals and print as integers when they round to an integer at that precision. Counts use integers. Positive values below the display precision are shown as an upper bound, for example \code{<0.0001}. Displayed F-statistics above 1000 are shown as \code{>1000}. Stored numeric values remain unchanged.}
\item{top,bottom}{Nonnegative integer numbers of leading and trailing rows to display
in each table. \code{NULL} uses \code{st.options()$print.top} and
\code{st.options()$print.bottom}, initially 10 and 5. Tables with no more than their sum are shown
in full. Otherwise an explicit omitted-row count separates the two parts.
Either argument may be \code{Inf} to show every row. Zero suppresses that part.
The model summary uses separate settings for the number of covariance blocks.}
\item{char.per.line}{Maximum output line width (integer 40--10000).
\code{NULL} uses \code{st.options()$print.char.perline}, initially 75.
The printer first reduces spacing, then decimal places (to at least two,
or one if explicitly requested), then abbreviates labels with a printed key.
Probabilities and existing numerical bounds keep the requested precision.
If needed, it splits wide tables into column blocks with repeated row numbers
and the first column. Extremely wide numbers use wrapped row records. Prose and
formulas also wrap. UTF-8 widths are measured conservatively in bytes without
splitting a character. Printing does not change stored numbers or labels.}
\item{effects}{\code{"fixed"} and \code{"random"} return the corresponding
cached test table. For \code{"random"}, a model without random terms returns
a zero-row table. Random terms can appear without a computed test.
Both \code{"both"} and \code{"all"} return one data frame of
stored fixed, random and MANOVA tests with the same columns for every fit method.
Its \code{scope} distinguishes individual traits, joint trait hypotheses and
whole terms. The reference distribution, degrees of freedom, test statistic,
\code{p} and denominator source are recorded where available. The
\code{statistic} and \code{calculation} columns retain the original MANOVA
measure and its F reference. The value \code{"manova"} returns only the
stored MANOVA table for one fitted model. Its columns are \code{term},
\code{method}, \code{statistic}, \code{DF1}, \code{DF2}, \code{F}, \code{p}
and \code{calculation}. The print method omits intercept tests by default
and identifies the term in the heading when a single term is displayed.
A fit without stored MANOVA tests gives an error. A second fitted
\code{st.lmm} object starts a multi-model comparison.}
\item{newdata}{Optional data frame for prediction. Fixed-factor levels and
retained transformations follow the fitted model. Only the predictors and
grouping variables used by the requested prediction are required. A random-only
prediction does not require fixed-effect predictors or offsets.}
\item{se.fit}{The prediction list always contains a \code{fit} matrix. TRUE
adds a \code{se.fit} matrix and, for ML and REML, a logical \code{estimable}
matrix. All three matrices have one row per observation and one column per
trait, including fits with one trait. For moments fits with standard errors,
the list also contains \code{DF} and \code{residual.scale}. A requested
moments interval appears separately in an \code{interval} matrix with
columns \code{lwr} and \code{upr}.
The common printer displays univariate diallel overlay predictions together in
one table. Diallel overlay results include readable cross labels such as
P2:P4; these are display metadata, and the returned numerical components
are unchanged.
For ML and REML fits, standard errors describe the prediction error of the
requested latent mean, including fixed-effect estimation when those
coefficients are estimated. They exclude future residual noise and
covariance-parameter uncertainty. Further prediction details appear below.}
\item{allow.new.levels}{TRUE allows previously unseen independent random-effect
levels with a predicted random contribution of zero. Their fitted prior
variance remains in conditional and random-only prediction uncertainty.
For a previously unseen covariance stratum, only point prediction with a
zero random contribution is available. Levels absent from a supplied known
relationship matrix remain errors; known parental relationships can define
new hybrid combinations.}
\item{new.levels}{Optional string convenience policy for unseen grouping
levels. \code{"error"} is equivalent to \code{allow.new.levels=FALSE};
\code{"population"} is equivalent to \code{allow.new.levels=TRUE} and gives
an unseen level zero random contribution. \code{NULL} preserves the existing
\code{allow.new.levels} argument and its default.}
\item{re.form}{For ML and REML conditional or random-only predictions,
NULL includes all fitted random coefficients. NA, \code{~0}, \code{~-1}
or \code{character(0)} includes none. A character vector selects complete
covariance blocks by their exact names in \code{st.lmm.components(object,
"random")}. A one-sided random formula selects exact fitted coefficients,
for example \code{~1|group} or \code{~0+x|group} from a fitted correlated
intercept and slope. Parenthesized sums such as
\code{~(1|group)+(1|block)} are also accepted. Selected coefficients must
have exactly the fitted basis and coding. The default marginal prediction
excludes all random effects; a nonempty explicit selector requires
\code{type="conditional"} or \code{type="random"}.}
\item{newparams}{Optional physical parameters for ML or REML prediction:
a named list with \code{fixed} and/or \code{covariance}. The latter is a
named list with \code{random} and/or \code{residual}, each naming every
covariance block in that section. Fixed coefficients must be a finite named
vector containing every fitted model-matrix column, including aliased
columns. Covariances must be named symmetric positive-semidefinite matrices
in the original coefficient units and must obey the fitted covariance
family. Omitted sections retain fitted values. Names, dimensions and family
constraints are checked, and coordinates are aligned by name. Further details below describe
the distinction between supplied and estimated fixed coefficients.}
\item{...}{Additional method arguments. For univariate fixed-effects fits
using \code{method="moments"}, \code{anova} accepts further nested fits.
For univariate ML or REML fits, \code{anova} accepts further fitted models
and \code{drop1} performs fixed-term deletion; further details appear in \code{\link{st.lmm.compare}}.
Fixed-effect subset comparisons use the \code{selection} argument of
\code{st.lmm}. Additional subset-summary arguments are described below.}
}
\details{
\subsection{Half diallel and parental overlay}{
\code{random=~overlay(mother, father)} creates one shared parental effect
vector. Each cross receives the sum of its two parental effects; a self
cross has loading two. The parental union is shared across both roles,
with one estimated GCA variance. An optional parent relationship matrix
is specified as \code{ucov=list(overlay=A)} or under the exact overlay term name.
Additional parent levels in this matrix are retained for prediction.
\code{\link{st.diallel.cross}} constructs an unordered cross factor
that can enter as a separate random term for SCA. Reciprocal and maternal
effects require separate specification.
The printed variance-component table reports the shared overlay variance
once, rather than expanding it into a repeated diagonal matrix. The
stored covariance matrix and fitted model are unchanged.


Univariate overlays can be standalone, interact with ordinary directly named
grouping factors, or be stratified by \code{at()}. The parental union is
common to both roles and all grouping levels or strata. A supplied parental
relationship matrix applies to every group and stratum; unused parent levels
in that matrix remain available for prediction. Parent roles may have different
factor level orders, and reversing the roles leaves the model unchanged.
A self cross retains loading two in every active stratum.

\tabular{ll}{
\strong{Random formula} \tab \strong{Covariance interpretation} \cr
\code{~overlay(mother,father)} \tab One shared parental-effect vector. \cr
\code{~overlay(mother,father):env} \tab Independent parental-effect vectors
in each environment, with one common variance. \cr
\code{~at(env):overlay(mother,father)} \tab Independent environment strata,
each with its own parental variance. \cr
\code{~at(env,"E2"):overlay(mother,father):year} \tab Parental effects only
in E2, independent across years, with one variance shared across years. \cr
\code{~overlay(mother,father):env:year} \tab Independent parental effects
for environment-by-year groups, with one common variance. \cr
}

Writing \eqn{h_i} for the two-parent incidence vector of observation \eqn{i},
a grouping interaction has covariance
\eqn{\sigma^2 h_i^T A h_j I(g_i=g_j)}. Its coefficient covariance is the
product of an independent grouping identity and the parental covariance.
An \code{at()} term adds the stratum indicators and estimates a separate
\eqn{\sigma_s^2} for each selected stratum. Multiple \code{at()} factors
form their selected Cartesian strata. Ambiguous combined stratum labels
are rejected and require distinct level labels before fitting. No correlation
across grouping levels is inferred. A common interaction variance alongside
all its separate stratum variances defines redundant covariance directions
and is rejected. Covariance-direction checks include fixed-effect projection
for REML. Active fixed/shared constraints use the existing reduced checks.

The common parental relationship is specified by \code{ucov=list(overlay=A)},
the exact overlay expression name, or its complete interaction term name. These
keys identify parental covariance; they do not specify covariance between
environments or years. A different matrix for each selected stratum is not
inferred. The representation stores a dense observation-by-parent basis and
parental covariance, so memory grows with observations and parent count.

\code{\link{st.predict}} with the complete fitted overlay block name returns
intercept plus parental GCA and its full joint prediction-error covariance;
\code{"0 + "} requests GCA alone. Grouped results list parents within each
independent group. Factor-level labels use byte-length prefixes to avoid
ambiguous concatenation, followed by the parent label. Stratum-specific blocks
are selected by their exact fitted names, available from
\code{names(st.lmm.components(fit, "covariance")$random)}.
Other fixed predictors are at their intercept reference; this is not an
average of fixed environment effects. Conditional \code{predict} with cross
data is the direct interface for complete cross means and their uncertainty.

Fitting, objective gradients, inference, refitting, updates, serialization,
training/future cross prediction and simulation use the same layout. New
independent grouping levels require \code{allow.new.levels=TRUE}; their
parental effects have zero conditional mean and the fitted prior covariance.
New parents must already be represented by the retained parental union or
known relationship matrix. New \code{at()} strata do not acquire new variance
parameters automatically.

Multivariate overlays, residual overlays, multiple distinct overlays in one
interaction, trait modifiers and units within overlay interactions, random
coefficient formulas and correlated coefficient products remain unsupported
and fail explicitly. Separate additive overlay terms are allowed when their
covariance directions are identifiable. No automatic maternal, reciprocal or
cross-environment correlation is introduced.
}
\subsection{Series of earlier fits}{
\code{series=list("env", fit1, fit2, fit3)} constructs the response data
from earlier fitted models without refitting them. The fixed genotype
factor is inferred as the single fixed predictor of the first fit shared
with the new fixed/random formulas, excluding the series factor. It must
be fixed and categorical in every earlier fit. Adjusted means and their
full model-based covariance use the compiled least-squares means machinery,
with equal nuisance-factor averaging and numeric covariates at its default
reference values. Nonestimable means are rejected. Transformed responses,
moment fits and model selection are outside this interface. Shared-predictor
ambiguity can be resolved by an explicit genotype in the reference declaration.

A multivariate second stage has a response formula such as
\code{cbind(y1,y2)~env+gen}. Each earlier fit must use those same untransformed
response names, with a fixed genotype factor. Trait order is matched by name.
All trait-specific genotype means and their full joint covariance are inherited;
cross-trait and cross-genotype covariance within an earlier fit is retained.
Covariance between independent earlier trials is zero. Missing first-stage
trait cells and incomplete genotype sets are allowed when every requested
within-trial genotype-by-trait mean is estimable. Otherwise extraction fails
explicitly. A row contributes to provenance and observation-ID checks once
when any trait is observed. The automatic mean reference grid uses the rows contributing any observed
trait, equal nuisance-factor averaging, and the usual numeric reference values.
Explicit reference grids and weights can replace that choice, as described
below. The reference definitions determine the targets when covariate distributions differ.

The second-stage data have one row per earlier-trial genotype mean and one
response column per trait. Known covariance is in trait-major order, with
all data rows for trait one followed by all rows for trait two. Its row and
column names identify trait cells as described below. Subsetting preserves
this order on both covariance axes. Refitting, updating, serialization,
training-data simulation and multivariate \code{predict} retain this covariance.
The existing univariate-only restriction in \code{st.predict} remains.
Ordinary model inference conditions on estimated first-stage covariance
parameters. For univariate independent trials, \code{st.lmm.bootstrap} with
an explicit \code{series} declaration can refit both stages; its current
scope and generation assumptions are documented there. Ordinary inference does
not account for first-stage covariance-parameter estimation uncertainty. Cross-fit dependence must be
absent or represented by the explicitly supplied joint covariance described below.

Additional trial-level predictors referenced in the new fixed or random
formulas are copied from each earlier fit's retained data and repeated for
its genotype means. Every predictor must be a nonmissing scalar column,
constant within each trial, with consistent categorical (factor/character)
or numeric interpretation across fits. Factors use their labels, so local
level codes need not agree. Numeric columns retain their numeric meaning;
Categorical years are represented by a factor. Predictors varying within a trial are
rejected rather than silently averaged.

\strong{Known cross-fit target covariance.}
The default treats errors in distinct earlier fits as independent.
\code{control$series.covariance} specifies the full known numeric covariance
matrix of all inherited target cells on the absolute response scale. It is not a correlation matrix and has no estimated multiplier.
Both axes must name every target cell exactly once; their orders may differ.
The \code{cells} entry of \code{st.lmm.components(object, "series")}
contains the complete cell declaration. Its columns are
\code{label}, \code{trial}, \code{genotype},
\code{trait}, and \code{data.row}. Cells are in trait-major order, before any
second-stage row subset. The labels encode the UTF-8 byte lengths and values
as \code{length(trait):trait:length(trial):trial:length(genotype):genotype}.
The returned labels identify cells without ambiguous separators.
A default fit can expose the declaration; its independence-based inference is
not valid evidence of independence when the caller knows dependence exists.

The supplied matrix must include all original target cells, even if the
second-stage call selects fewer rows. Every within-fit block, including
cross-trait and cross-genotype entries, must agree with the extracted covariance
to absolute tolerance \eqn{10^{-8}} on the earlier-target correlation scale.
Those blocks retain the exact extracted values. The full matrix must be finite,
symmetric and positive semidefinite, checked after normalization by the earlier
target standard deviations. Earlier target variances must be positive and finite.
Observed marginal covariance still must be positive definite; no jitter or
ridge is added. Tiny accepted cross-block asymmetries are averaged. Reference
grids and weights define the targets before this alignment and validation.

This option handles externally known dependence, for example a shared
calibration error across distinct observations. It does not infer or estimate
cross-fit covariance from matching IDs or propagate uncertainty in estimated
first-stage covariance components. Correct cross-fit values and meaningful
target definitions remain the caller's responsibility. Copied fits, overlapping
source provenance and duplicate observation IDs remain rejected: a supplied
matrix does not authorize reusing observations in separate earlier fits.
The \code{covariance.provenance} entry of the series component records the
supplied-covariance source, within-fit tolerance and pre-subset scope. The
cell mapping and source are retained through refitting, updating,
serialization and simulation. The \code{rcov} entry retains its fit-time
dimnames, when present. The \code{series.covariance} control uses the
complete set of labels in the \code{cells} entry.

\strong{Explicit target/reference grids.}
An explicit reference is specified by
\code{control=list(series.reference=list(genotype="gen", grid=grid,
weights=weights))}. The \code{genotype} name is optional when the usual
single shared fixed-factor inference is unambiguous; declaring it resolves
ambiguity when additional shared covariates occur in both stages. It must
name a fixed categorical factor in every earlier fit and occur in the
next-stage fixed or random formula. The multivariate virtual \code{trait}
modifier cannot be the genotype target.

\code{grid} is either one common data frame for all trials, or a complete
named list of data frames indexed by the resolved trial labels. Each grid
contains the genotype column and every raw predictor used by the earlier
fixed formula, excluding the multivariate virtual trait. No extra columns
are accepted. Columns are scalar, nonmissing and finite where numeric;
genotype targets use factor or character labels. All values are interpreted
through the earlier fit's retained design schema, including learned
transformations and factor contrasts. A common grid with \code{x=20}, for
example, requests the same covariate reference in all trials rather than
their different observed means. Trial-specific grids can accommodate
trial-specific nuisance-factor labels or explicitly different references.

The requested genotype targets and their order are the distinct genotype
labels in first-appearance order in each grid. The grid is used exactly as
supplied; no Cartesian expansion or automatic nuisance-factor averaging is
added. Replicate rows can express the desired joint distribution of nuisance
factors and covariates. This permits valid nested targets that the automatic
grid rejects. Unknown design levels and nonestimable weighted linear targets
fail explicitly. Individual rows need not separately be estimable when their
weighted linear combination is estimable; the caller must declare meaningful
reference combinations.

\code{weights=NULL} gives equal weight to grid rows within each genotype.
Otherwise one finite nonnegative numeric weight is required per row. A named
vector is aligned independently to the grid's row names; names must match
exactly once. A complete named list of numeric vectors or \code{NULL} values
can supply trial-specific weights. Weights are normalized within genotype,
not across genotype means or trials; every target needs positive total
weight. Uniform rescaling within a target leaves the result unchanged.
Zero-weight rows still undergo grid and design validation. The full joint
covariance is computed from the same weighted linear functions as the means.
These are reference-distribution weights, not residual precision weights.
A trial's grid and weights apply to all response traits; trait-specific
reference distributions are not part of this interface.

For example, \code{expand.grid(gen=c("G3","G1"), block=c("B1","B2"))}
with a declared common \code{x} column and row weights describes explicit
weighted nuisance-factor averages for the two targets. All fixed predictors required by each earlier model must be present. Per-trial grids need not request
the same genotype set, but the second-stage model still needs sufficient
connectivity and estimability.

When an additional next-stage predictor is present in the reference grid,
its value must be constant over positive-weight rows within each genotype
target. That target value is copied into the second-stage data. A varying
nuisance factor cannot silently be represented by one category or averaged
numeric value in a new formula. Other additional next-stage predictors retain
the existing requirement of being constant within the earlier trial.
The \code{reference} entry of the series component retains the declared grids
and weights. Its \code{genotype} entry records the resolved genotype column.
These declarations precede any second-stage row subset. Refits,
updates and response simulations retain them without re-extracting means or
refitting the earlier models. With no explicit reference declaration, the
automatic-grid behavior remains in use.

For repeated locations across years, each independent location-year trial
has a separate fit, common loc and year identifiers, and a unique env trial label. A second-stage fixed formula y~env and random formula
~gen+gen:loc+gen:year+gen:env separates main genotype, persistent
genotype-location, genotype-year and transient genotype-trial effects.
Incomplete grids are allowed, but sufficient connectivity is needed to
identify the requested effects. An independent genotype-trial variance can
be redundant with an additional units variance, as checked below. The same
covariance check applies to an equivalent gen:loc:year interaction.
A supplied genomic covariance for gen does not automatically apply to its
interaction terms. No temporal correlation or future-year forecasting is
created merely by retaining year.

Each earlier fit must represent distinct observations from a trial. Errors
are assumed independent across fits unless their known joint target covariance
is supplied explicitly. No cross-fit covariances are inferred. Named fit entries supply environment labels when
the environment column is absent. If that column is retained, it must have
one nonmissing value on the observations used by the earlier fit, and a
supplied entry name must match it. Named entries cannot combine multiple
environments into one trial. Duplicate environment labels are rejected.
The response, omission indicators, fixed-design dimensions and retained
row map are checked before adjusted means are extracted. An unused matrix
column does not determine the observation count.

New top-level fits carry an opaque provenance token which is preserved by
serialization, \code{st.lmm.refit}, and \code{update} when the original data
are retained. Their original row maps detect shared contributing observations.
Copies and overlapping refits cannot be counted as independent trials;
disjoint subsets of a common source can be supplied as separate trials.
Separately fitted independent
trials are not equated merely because their values happen to agree.
Identical legacy snapshots without a provenance token are rejected when
independence cannot be checked. Independence can be checked after refitting independent trials with the
current version or through genuinely disjoint observation IDs as follows.

With \code{control=list(series.id="plot.id")}, every earlier fit must retain
the scalar column \code{plot.id}. Its IDs must be nonmissing, nonempty and
unique among the observations used for that fit. IDs shared by different
fits are rejected, including after serialization. Factors are compared by
labels; numeric IDs use a full-precision representation. Rows with
missing responses or omitted model inputs do not contribute observations
to this check. Relabelling shared observations does not establish
independence; shared evidence requires a joint analysis.

The \code{independence} entry of the series component records the trial labels,
observation counts, source provenance, ID column and whether ID disjointness
was checked or independence remains the caller's assumption. Disjoint IDs
do not prove that trial errors are independent: shared calibration and
other cross-fit dependence require an externally known joint target covariance
or a joint analysis. The response
and its covariance are assembled in matching order, and existing subset
and covariance-alignment handling applies thereafter. The data and rcov
arguments cannot accompany series.

The printed residual description explicitly identifies the supplied covariance
and absence of an additional units variance when cov is NULL. The inherited
(or supplied) covariance is reported separately from estimated components,
with its scale fixed at one; this is not a common error variance of one.

The inherited covariance is the block diagonal matrix of earlier adjusted
mean covariances, with coefficient fixed at one. It conditions on the
estimated earlier covariance components. Ordinary units covariance adds
an estimated independent GxE variance. For scalar variance components,
identifiability is checked using the actual covariance matrices on the
retained observations. REML checks project out the fixed-effect design.
If the omitted default units component is the sole redundant covariance
direction, it is omitted automatically. Explicit \code{cov=~units} is
rejected when it is redundant. Identifiable alternatives include
\code{cov=NULL}, omission of a redundant term, or fixed/shared constraints.
A random environment main effect alone does not suppress units. A genomic
GxE kernel can coexist with units when their projected covariance directions
are independent. Raw matrix rank before accounting for fixed effects is not
sufficient: a compound-symmetric GxE kernel can become proportional to units
after fixed environment effects are removed. Fixed/shared constraints are
checked in their independent feasible scalar covariance directions.

For the shared genetic random term, \code{st.predict} on a series fit adds
the equally averaged fixed environment mean instead of the reference
intercept. Its SE uses joint prediction-error covariance. The returned
prediction has a \code{reliability} attribute, also shown in its summary,
equal to genetic BLUP variance divided by prior genetic variance. This
reliability excludes uncertainty in the added environment mean. Genotyped
levels represented in the supplied relationship matrix are retained for
prediction, including those without phenotypic means. Missing relationship
matches remain errors under the ordinary ucov contract; no lines are
silently discarded.

The maintained v-lmm companion contains the three-environment example.
This implements the homogeneous GxE
stagewise model of Endelman (2023), doi:10.1007/s00122-023-04298-x, using
full earlier-stage covariance propagation as in Piepho et al. (2012),
doi:10.1002/bimj.201100219. It does not implement StageWise's spatial
splines or factor-analytic location covariance.
}
\subsection{Known residual covariance}{
With \code{rcov=C1,cov=NULL}, the residual covariance is exactly \code{C1},
without estimating a multiplier. This propagates the complete covariance
of first-stage adjusted means into a second-stage model.
\code{as.matrix(vcov(means))} from \code{st.lsm} provides the covariance;
its row names identify the corresponding rows of the second-stage data.
Diagonal precision weights discard these between-mean covariances.

The matrix must be finite, symmetric and positive semidefinite; the
observed marginal model covariance must be positive definite. No automatic
ridge or jitter is added. If either axis is named, both axes must have
unique names matching \code{rownames(data)}. They are aligned independently
before applying \code{subset} and missing-value exclusions to both axes.
Without names the matrix is positional in the original data row order.
For \code{cbind(y1,y2,...)} with \eqn{n} input rows and \eqn{p} traits,
\code{rcov} is an \eqn{np\times np} matrix on the original trait scales.
Its positional order is all \eqn{n} rows for trait one, then all rows for trait
two, and so on. Named axes use explicit trait-cell labels:
\code{paste0(nchar(trait, type="bytes"), ":", trait, ":", rowname)}.
The byte-length prefix makes labels unambiguous even when names contain
colons. For example, trait \code{y1} on row \code{plot7} is
\code{"2:y1:plot7"}. Both axes must contain the complete set exactly once;
row and column permutations are aligned independently. Names use UTF-8
trait labels. Trait names must be nonempty, nonmissing and distinct.
Missing responses select the observed principal submatrix without replacing
unobserved responses. A row subset selects its entries within every trait,
including every cross-trait covariance block. The covariance must include
all original rows and traits before this selection.

An additional \code{cov} formula adds its estimated covariance components
to the supplied matrix. The fit conditions on the supplied covariance;
uncertainty in its first-stage estimation is not propagated.

The matrix is retained for updating, refitting and training-data
simulation. Future residual simulation requires future covariance and is
therefore rejected by the current newdata interface. Latent-mean newdata
prediction remains supported. Known \code{rcov} is not supported by the
method of moments or with precision weights.
}
\subsection{Row selection and missing values}{
Response and formula expressions, including data-dependent transformations,
are evaluated on the original data before applying \code{subset} and the
explicit missing-value policy. This follows \code{stats::model.frame}; a
polynomial or spline is therefore calibrated on the full input even when
only a subset is fitted. An expression that cannot be evaluated on the full
input must first be supplied as a precomputed data column.

Explicit \code{na.omit} and \code{na.exclude} remove a common row whenever
any referenced response, predictor, grouping variable, offset or precision
weight is missing. For a matrix response, a missing trait removes that row
for every trait. Unreferenced data columns do not affect this policy.
\code{na.fail} rejects missing model inputs after subsetting.
\code{na.pass} retains missing responses as prediction rows but requires
finite complete predictors, offsets and weights. Nonmissing infinite
responses are invalid. With an explicit policy, retained predictors and
offsets must be finite. With the default NULL, missing responses and the
historical predictor-omission behavior remain available.

For mixed models and structured residual models, every retained precision
weight must be finite and strictly positive after the explicit NA policy.
Ordinary weighted regression retains its existing exclusion of zero and
non-finite weights; negative finite weights remain errors. Subset comparison
additionally uses its common complete analysis sample.

Internal matrices and \code{model.frame} contain retained rows. With
\code{na.exclude}, public \code{fitted}, \code{residuals} and predictions
without \code{newdata} insert NA values at removed positions in the selected
row order, including repeated selected rows. A new-data prediction uses only
its new-data rows. The model retains original row IDs, unique selected row
names, retained positions and separate NA, weight and common-sample omission
maps in \code{model$input}. Further details appear in \code{\link{st.lmm.model}}.
}
\subsection{Comparing fitted models}{
\code{anova(m1, m2, ...)} compares fitted univariate likelihood models.
The default \code{refit=TRUE} converts REML fits to ML and reports this
convention in the output. \code{drop1(m)} deletes eligible fixed formula
terms while preserving hierarchy. Observation compatibility, nesting and
variance-boundary qualifications are described in
\code{\link{st.lmm.compare}}.
}

\subsection{Offsets, contrasts and retained model information}{
An additive \code{offset(expression)} in a univariate fixed formula supplies
a known contribution to the response mean. Multiple offsets are added;
they do not create estimated coefficients or model-matrix columns. Estimation
uses the offset-adjusted response, while fitted values and predictions are
returned on the response scale. New-data predictions evaluate the retained
offset recipe in the new data. Offsets cannot be random-coefficient terms or
parts of formula interactions. Least-squares means include the known offset evaluated
on their reference grid, and trends include its finite-difference contribution.

The explicit \code{offset} argument and formula offsets are stored separately
and each contributes once. A scalar argument offset, or an expression based
on data columns, can be reused for new-data predictions. Learned scaling and
polynomial or spline recipes retain their training calibration. A training-only
external offset vector cannot be safely matched to new data, even if the row
counts happen to agree; an explicit prediction \code{offset} is required in
that case. Least-squares means and trends also evaluate scalar or data-derived argument
offsets on their reference grid. A training-only vector without a grid recipe
is rejected by \code{st.lsm}; no arbitrary reference offset is assumed.

Factor contrast matrices are stored with their level mappings. Later changes
to factor attributes or global contrast options do not change the fitted
coding used for prediction. No-intercept terms use full indicator coding
where required by the formula. Reduced or redundant user contrast matrices
specify their actual model space; fixed-effect aliases are reported by the
likelihood and fixed-only moments paths. For reduced or redundant coding,
Type III construction centers the supplied columns over factor levels.

\code{model.frame}, \code{terms} and \code{weights} expose retained model
information. \code{update} constructs a new fit using retained data and
settings; \code{update(object, evaluate=FALSE)} returns its call. Further details appear in
\code{\link{st.lmm.model}} for evaluation and persistence details.
}
\subsection{Precision weights}{
For a supported structured residual covariance \eqn{R_0}, the supplied
precisions define \eqn{R = W^{-1/2} R_0 W^{-1/2}}, where
\eqn{W=\mathop{\mathrm{diag}}(w_i)}. Random-effect covariance is unchanged.
The transformation also applies with known random-level covariance matrices,
bitpacked SelectionTools marker data sets and SelectionTools hybrid data sets.

When all residual covariance parameters are estimated, multiplying every
weight by the same positive constant rescales those parameters by that
constant and leaves the fitted observation covariance, coefficients,
predictions and likelihood unchanged. The reported residual variance therefore
refers to a weight of one. Missing-response prediction rows in mixed models
remain in the data and require valid positive precisions.
}
\subsection{Fixed-effect subset comparison}{
Supplying \code{selection} compares subsets of the fixed regression part of
a univariate model. The selection unit is explicitly one original column of
\code{data}, not one column of the expanded model matrix. All contrast columns
of a categorical factor enter or leave together. Transformations using the
same source column also enter or leave together; for example, \code{x} and
\code{I(x^2)} form one selection unit when both occur in the full formula.
An interaction or other term involving several source columns is retained
only when all of those source columns are selected. Selection does not add
terms absent from the supplied full formula. This grouped search deliberately
differs from procedures that select individual factor contrasts.
The full fixed formula must reference at least one candidate data column.
At most 63 original source columns are supported by the subset search.

The fitted object continues to represent the full model specified by
\code{fixed}, \code{random} and \code{cov}. Subset results are attached to
\code{object$selection}; no candidate automatically replaces the fitted model.
The default \code{summary} shows the comparison when a search is present.
\code{summary(object, selection=FALSE)} provides the ordinary full-model
summary. The user may inspect the comparisons and then make another
\code{st.lmm} call with the chosen formula. The \code{id} arguments provide
access to retained candidates without changing the original fitted object.
Refits for candidate extraction use the requested estimation method and are
performed when needed, rather than fitting every retained candidate for full
inference during the search.

For \code{method="moments"}, selection requires fixed-only ordinary or
weighted regression with \code{cov=~units}. For \code{method="ML"} or
\code{method="REML"}, supported univariate random effects and residual
covariance structures may remain in the model. Only fixed regression columns
are selected. The random and residual structures, supplied covariance
matrices and observation set remain unchanged across candidates, while
unknown covariance parameters are re-estimated for each candidate.

For likelihood-based fits, fixed-effect subsets are compared using ML
likelihoods. In particular, a
call with \code{method="REML"} fits the specified full model by REML but
switches to ML for subset comparison. Both methods are identified in the
printed output and stored result. REML likelihoods are not used to rank
models with different fixed-effects designs. A subsequent ordinary
\code{st.lmm} call can fit the user's chosen formula by REML.

The following named components are accepted in \code{selection}:
\describe{
\item{\code{search}}{One of \code{"exhaustive"} (default),
\code{"forward"}, \code{"backward"} or \code{"seqrep"}.
Exhaustive search examines all allowed subsets. Forward selection adds a
source column at each step; backward selection removes one. Sequential
replacement improves visited subsets by exchanging source columns.}
\item{\code{nvmax}}{Largest number of selected original data columns;
default 8, limited by the available columns. The intercept does not count
towards this limit. Zero requests the model without any candidate source
column, retaining the intercept when it is part of the fixed formula.
Ordinarily the smallest reported model contains one more selected column
than the forced-in set. If \code{nvmax} equals the number of forced-in
columns, the model containing exactly that set is retained instead; this
also permits \code{nvmax=0} with an empty forced-in set.}
\item{\code{nbest}}{Integer from 1 to 10000 giving the number of best candidates retained for
each combination of source-column count and estimable fixed-effect rank;
default 1. Retaining separate ranks prevents a wider factor from displacing
a candidate preferred by a penalized criterion. A comparison summary can
therefore contain more than \code{nbest} rows for a given source-column count.}
\item{\code{nested}}{Logical; defaults to \code{nbest==1}. For forward and
backward search, TRUE reports the single nested path, using one retained
path model at each size regardless of \code{nbest}. With FALSE, forward
search also considers visited additions and initial source-column prefix
models; backward search considers initial and changing path-prefix models.}
\item{\code{force.in}, \code{force.out}}{Source-column names, indices in
\code{object$selection$groups}, or logical vectors with one entry per group,
identifying columns that must be included or excluded. They cannot overlap.
These controls refer to complete data columns, not individual factor levels
or contrast columns. Intercept presence is controlled by the fixed formula.}
\item{\code{really.big}}{Logical; default FALSE. Explicitly permits an exhaustive search
whose estimated candidate count exceeds the large-search threshold:
100000 for ordinary regression or 10000 for a fitted covariance model.
It raises the default \code{max.models} limit but does not remove that limit.}
\item{\code{max.models}}{Maximum number of candidate evaluations; default
100000, or 10000000 when \code{really.big=TRUE}. Ordinary-regression
reference evaluations count towards this limit: one for the full model,
or two when \code{force.out} distinguishes the allowed and requested full
models. Intermediate reduced-model evaluations used to score sequential
replacements also count towards the limit.
Heuristic searches cache at most 1000000 distinct candidates without
reducing the requested evaluation limit. An exhaustive request exceeding
the evaluation limit is rejected before the search; a heuristic reaching its
limit is reported as incomplete.}
\item{\code{criterion}}{Optional criterion preference for the comparison
summary, supplied as a character vector. NULL uses all applicable criteria.
All applicable fit statistics are calculated even when a preference is given.
This preference controls the summary view and criterion winners; it does not
change the greedy search path, select or replace the final fitted model.}
}

Ordinary regression candidates are retained by residual sum of squares within
each source-column count and fixed-rank class. Candidates with fitted covariance
structures are retained by ML log-likelihood within these classes. Forward, backward and sequential
replacement searches are heuristic and do not guarantee the globally best
subsets. Their paths are not stopped merely because a penalized criterion
worsens. Exhaustive search is complete only if all required candidates were
evaluated successfully; fit failures and an incomplete search are reported.
Covariance optimization convergence is checked separately from coverage of
the subset search.
The retained candidate table is limited to 5000000 slots over combinations
of source-column count, fixed-effect rank and \code{nbest}; requests exceeding
that implementation limit must reduce \code{nvmax} or \code{nbest}.

Subset summaries accept the following additional arguments through
\code{...}:
\describe{
\item{\code{criterion}}{NULL uses \code{selection$criterion} when supplied,
otherwise all applicable criteria. A character vector selects one or several
criteria from \code{"rss"}, \code{"rsq"}, \code{"adjr2"},
\code{"Cp"}, \code{"BICm"}, \code{"AIC"}, \code{"BIC"} and
\code{"logLik"}. Only the last three apply to fitted covariance models.
The aliases \code{"RSS"}, \code{"r2"}, \code{"R2"}, \code{"cp"},
\code{"aic"} and \code{"bic"} are also accepted. The summary contains
the union of the best requested candidates per source-column count for each
criterion. Consequently, different criteria can contribute different rows
at the same size. A requested criterion that is undefined for every retained
model produces an error; the automatic criterion set leaves such a winner
unavailable.}
\item{\code{all.best}}{Logical; TRUE (default) uses the search's
\code{nbest} for each requested criterion and size. FALSE uses one best model
per criterion and size. It does not change the stored search results.}
\item{\code{matrix}}{Logical; TRUE (default) constructs the summary's
inclusion table \code{outmat}. FALSE leaves that component NULL.}
\item{\code{matrix.logical}}{Logical; FALSE (default) represents inclusion
in \code{outmat} with stars and blanks. TRUE constructs a logical matrix.
The \code{which} component is always logical and includes the intercept
column when fitted; \code{outmat} omits that intercept column.}
\item{\code{DF}}{Optional positive scalar for ordinary regression. It
replaces \eqn{n-I} in the usual summary formulas, where \eqn{I} is 1
when an intercept is fitted and 0 otherwise. The corresponding effective
sample count is \eqn{DF+I}. This changes adjusted \eqn{R^2}, \eqn{C_p} and
\code{BICm}; it does not refit any candidate or change the actual-sample
ML log-likelihood, AIC or BIC. NULL uses the observed sample count.}
}

All candidates use one common analysis sample determined from the complete
specified model and weights. Missing or non-finite responses, predictors,
grouping values or weights are excluded consistently, and zero-weight rows
are excluded. Removing a candidate predictor does not add rows back to the
sample. A separate refit supplied with the original data can have a different
sample if its reduced formula no longer refers to a missing predictor.
Reproduction of a comparison score requires the same analysis data. The prepared
analysis data are retained in \code{object$data}; \code{selection$rows}
records the corresponding original row indices.

For ordinary regression, the comparison includes RSS, \eqn{R^2}, adjusted
\eqn{R^2}, Mallows' \eqn{C_p}, \code{BICm}, AIC, BIC and the ML
log-likelihood. Here \eqn{n} is the common number of observations,
\eqn{r} is the estimable fixed-effect rank including any intercept,
\eqn{RSS} is the residual sum of squares and \eqn{TSS} is the corresponding
total sum of squares. \code{BICm} is the normalized regression statistic
\deqn{BICm = n\log(RSS/TSS) + r\log(n).}
This reproduces the ordinary-regression statistic called \code{bic} by
\code{regsubsets}, without using that package. \code{BIC} retains the usual
likelihood definition. RSS uses the supplied precision weights.
For comparisons with an intercept, TSS is the weighted sum of squared
deviations from the weighted response mean. Without an intercept, the
comparison uses \eqn{TSS=\sum_i y_i^2}, without weighting, even when RSS is
weighted. This follows the legacy \code{regsubsets} normalization for
\code{rsq}, \code{adjr2} and \code{BICm}; the sum includes only the common
analysis sample, after excluding missing, non-finite and zero-weight rows.
This reporting convention does not change the weighted fits, RSS,
log-likelihood, AIC, BIC or the full fitted model's own weighted
\eqn{R^2} calculation. Undefined statistics,
including \eqn{C_p} when its reference variance is unavailable or nonpositive,
are reported as unavailable. Fitted covariance models report ML log-likelihood,
AIC and BIC; ordinary-regression \eqn{C_p}, adjusted \eqn{R^2} and
\code{BICm} are not assigned a mixed-model interpretation.

The \eqn{C_p} reference follows the \code{regsubsets} compatibility convention,
including when \code{force.out} is supplied. The quantity
\eqn{RSS_{\mathrm{requested}}} is the residual sum of squares of the full
specified fixed formula, before applying \code{force.out}, and
\eqn{r_{\mathrm{allowed}}} is the estimable fixed-effect rank of the model
containing every allowed source column after applying \code{force.out}.
Then
\deqn{s^2_{C_p} = RSS_{\mathrm{requested}}/(n-r_{\mathrm{allowed}}),
\qquad C_p = RSS/s^2_{C_p} - n + 2r.}
Thus the numerator and residual degrees of freedom can come from different
models when columns are forced out. This is a compatibility convention; with
no forced-out columns, both refer to the same full model. The reference models
are not limited by \code{nvmax}. For \code{summary(DF=...)}, the effective
sample count \code{DF+I} replaces \eqn{n} in both formulas.

The stored components \code{full.allowed.rss} and \code{full.allowed.rank}
describe the allowed full model; \code{full.rss} and \code{full.rank} retain
the same meanings. The components \code{cp.reference.rss},
\code{cp.reference.rank} and \code{cp.reference.DF} record the requested-full
RSS, allowed-full rank and residual degrees of freedom actually used for
\eqn{C_p}. The \code{full.allowed.*} and \code{cp.reference.*} components
are also available in the comparison summary, where \code{cp.reference.DF}
reflects any summary \code{DF} override.

The comparison stores inclusion indicators in \code{which}, source-column
names in \code{groups}, and retained formulas in \code{formulas}. Its
\code{statistics} data frame contains \code{id}, \code{size}, \code{rank},
\code{parameters}, \code{rss}, \code{rsq}, \code{adjr2}, \code{cp},
\code{BICm}, \code{AIC}, \code{BIC} and \code{logLik}. Here \code{size}
counts selected data columns, \code{rank} counts estimable fixed coefficients,
and \code{parameters} includes estimated covariance parameters for the
likelihood criteria. Summary objects also provide the statistic vectors
directly and a \code{winners} data frame with the criterion, model identifier
and value of each preferred retained model. Model identifiers denote retained
rows, not model sizes. Plotting displays an inclusion table ordered
by the requested criterion. A criterion can identify different preferred
models without rerunning the search because rank-specific candidates are
retained. Cross-validation, folds and resampling remain the responsibility
of the calling R code; the package supplies fitting and prediction.

Search diagnostics include \code{evaluated}, \code{failed}, \code{complete}
and \code{optimal}. \code{complete} means that the chosen search procedure
finished without failed candidate fits or a reached evaluation limit;
\code{optimal} additionally requires exhaustive search. The
\code{failures} component stores inclusion patterns and messages for the
first 100 failed fits. The full count remains in \code{failed}.
}
\subsection{Moment estimators}{
With \code{method="moments"}, fixed-effects models use least squares;
multivariate fixed-only models with one common fixed design also provide
MANOVA tests. Mixed designs require a balanced orthogonal layout.
Predefined designs retain their original default calculations. Other
balanced orthogonal crossed or nested categorical designs use the
general Cornfield--Tukey/Underwood path. Sources without an available exact or quasi-F test
are retained with status \code{"no exact test"} and unavailable F and p
values. Quasi-F tests are enabled by default in the general
calculations and use sums of two independent mean squares on
each side and Satterthwaite DF for both sides. Unsupported nonorthogonal
layouts and unresolved covariance inference raise errors; no likelihood fit is substituted. Linear functions combining
independent error strata use Satterthwaite DF by default; this does not
change the design-specific exact ANOVA F tests.
Negative variance-component estimates are retained without truncation.
Certified balanced categorical layouts use centred grouped projections for
source sums of squares and sums of products, with analytic EMS coefficients.
Integer cell-count checks verify the source intersections and orthogonality;
equal marginal factor counts alone are insufficient. This calculation keeps
the same moment equations, restricted interactions and inference methods.
The numerical orthogonal-basis calculation remains available for layouts or
rank tolerances outside that certificate. Mixed categorical fits and
fixed-only fits using \code{control$design="general"} record the selected
calculation in \code{object$moments$calculation}.
The stored estimation details are in \code{object$moments}, including
residual sums of squares, mean squares and, where applicable, variance
components and expected mean squares (EMS). Printing shows EMS in a separate
section immediately after the residual section. Displayed formulas use
\code{s2} and labels such as \code{EMS(Intercept)}. Full factor names are
retained when the line fits; unique short labels and an abbreviation key
are used only for expressions exceeding \code{char.per.line}. EMS equals
signs are aligned. EMS component columns are also aligned
when the complete layout fits; otherwise only the equals signs are aligned.
Test columns are displayed in the order \code{DF1}, \code{DF2}, \code{SS}, \code{MS},
\code{F}, \code{p}, where present. The
numerical test and residual tables retain their full stored \code{EMS}
expressions. Directly extracted test tables, MANOVA tables and numeric
estimation components use the common compiled printing controls, including
\code{digits}, \code{top}, \code{bottom} and \code{char.per.line}.
Moment fits support \code{control$rank.tol} and
\code{control$definite.tol}, plus the categorical analysis options below.
Likelihood, AIC and BIC are unavailable.

\describe{
\item{\code{design="auto"}}{Predefined paths are retained, with automatic fallback
for other balanced orthogonal designs. \code{design="general"} selects the
general path for univariate fits and multivariate mixed fits. Multivariate
mixed fits use common complete observations and a common fixed design for
the traits.
An intercept and positive residual DF are required. Factors must have complete
observations and balanced observed groups. Fixed nested factors, random
regressions, relationship matrices and nonorthogonal sources are unsupported.}
\item{\code{quasi.f=TRUE}}{A quasi-F search is performed by default when no
single exact denominator exists. The EMS under the null must match on both
sides. The actual source DF are stored in \code{sourceDF}, the approximate
numerator DF in \code{DF1}, denominator DF in \code{DF2}, and mean-square
combinations in \code{numerator} and \code{denominator}. P-values use full
precision DF, unlike GAD 2.0 which rounds these DF before evaluation.
Unresolved tests remain explicitly unavailable. Exact tests take precedence.
\code{FALSE} requests exact tests only.}
\item{\code{pool=character(), pool.method="pool"}}{Explicit sequential source
names; no automatic selection by p-value. \code{pool.method} is one method
or one per source. \code{"pool"} merges the selected source's SS and DF with
its single exact denominator, named \code{Pooled1}, etc.; \code{"eliminate"}
merges them with the residual. Each step removes that source component from
the EMS and reconstructs tests. Pooling requires a single exact denominator; elimination accepts an
available exact or quasi-F test. Untestable sources cannot be removed. Names in later steps refer to the current table.
\code{object$moments$pooled} contains the ANOVA, history, components and residual.
Multivariate mixed fits apply the same source sequence to all traits and
combine the full sums-of-squares-and-products matrices and their DF.
The pooled object also contains \code{SSCP}, \code{mean.products}, and
\code{strata}. The pooled source matrices are available as
\code{st.lmm.components(object, "sscp")$pooled}.
\code{anova} and printing use the pooled hypothesis table.
The fitted model,
\code{coef}, \code{vcov}, \code{fitted}, \code{residuals}
and \code{st.lsm} retain the original unpooled model and strata, stored in
\code{object$moments$anova}. Thus pooling does not silently change the
hypothesis represented by a subsequent coefficient contrast. Elimination uses
the updated residual DF in p-values (GAD 2.0 sometimes uses the original DF).}
\item{\code{components=FALSE}}{With TRUE, a components-of-variation table
for fixed and random sources is stored in \code{object$moments$components}. A tested source's
estimate is its test numerator minus denominator mean square, divided by its
own EMS coefficient. A source without a test has an unavailable summary
estimate, even when the full random-component moment equations are solvable.
Random variance estimates from those full equations remain separately stored
in \code{variance.components}. Negative estimates are retained, their SD is
unavailable, and their percentage is zero. Percentages follow GAD: shares of
the sum of nonnegative component SDs, not shares of variance. Multivariate
mixed fits report these summaries and percentages separately by trait.
Pooled summaries are stored in \code{object$moments$pooled$components}.}
\item{\code{cochran=FALSE}}{With TRUE, Cochran's C, its bounded p-value,
observed cell count, replication, largest-variance cell, cell variances and
means are stored in \code{object$moments$cochran}. Cells are grouped by
the factors in the model. Equal replication of at least two observations
per cell is required. Multivariate mixed fits use the same grouping for
every trait and return a named list with one diagnostic per trait.
Each p-value uses the GAD/Underwood F-based approximation, capped at one.}
}
\code{quasi.f=TRUE} and explicit pooling select the general path.
Component summaries and Cochran's C use the selected fitting path.
Mean comparisons continue through \code{\link{st.lsm}}; SNK is not implemented.

In a nested-model comparison, differences and test statistics for the
first (reference) model are displayed as blank cells. Their stored values
remain numeric \code{NA}. The exact F tests use the residual mean square
and denominator DF of the largest model in the comparison. The table
reports these denominator DF in \code{DF2}. With \code{test="none"},
the comparison reports changes in rank and residual sums of squares
without F tests.

Replicated regression uses the ordinary fixed formula and sequential
tests: \code{st.lmm(y ~ x + factor(x), data=d, method="moments", type=1)}.
The numeric term must precede the factor term. In \code{anova(m)}, the
\code{x} row tests the linear trend and the \code{factor(x)} row tests
lack of fit; the residual is pure error. With \eqn{g} distinct predictor
values and \eqn{n} observations, these have 1, \eqn{g-2} and \eqn{n-g}
DF, respectively. Both F tests use the pure-error mean square. Replication
and positive pure-error variance are required for inference.

Fixed-only moment fits use the design rank for residual and hypothesis DF.
Redundant columns, including those in \code{x + factor(x)}, are retained
in the coefficient coordinates and handled by a pivoted QR basis. Linear
functions are checked for estimability. Aliased coefficients are reported
as \code{NA}; a term with no estimable hypothesis has zero numerator DF
and unavailable F and p values. Predefined mixed moment designs still
require a full-rank fixed design.
Type I tests include the sequential intercept contribution before all
predictor terms. Its sum of squares is \eqn{n\bar y^2}; its F test uses
the fitted model's residual mean square, as do the later sequential terms.
This tests the overall mean and differs from an adjusted intercept test.
When every raw numeric coefficient is non-estimable, the model printout
omits the regression coefficient table; the stored \code{regression}
component retains the coefficient estimability results.

Multivariate mixed moment fits use common complete observations and the
same categorical fixed design for every trait. Ordinary categorical random
intercepts and \code{cov=~units} are available for balanced orthogonal
designs. The general path constructs a common set of strata from the
factor structure. Each trait's mean squares determine its exact or
quasi-F test. One example is the multi-environment layout
\code{random=~gen+env+gen:env+rep:env}, with genotypes crossed with
environments, equally replicated blocks within each environment, and one
observation per genotype/environment/replicate cell. Labels of replicates
are common across environments. With \code{fixed=cbind(y1, y2) ~ gen}
and \code{random=~env+gen:env+rep:env}, the same design estimates
least-squares means for each genotype and trait. \code{st.lsm(fit, ~gen)}
returns their estimates, SE, DF and confidence limits. With
\code{random.test="F"}, the general path uses exact denominators where
available and searches for quasi-F tests for the remaining terms.

Each orthogonal source gives a sums-of-squares-and-products (SSCP) matrix.
Dividing by source DF gives a mean-product matrix. The engine solves the
expected mean-product equations for every variance and covariance entry
using one common coefficient matrix, without likelihood optimization.
For the multi-environment layout with random genotypes,
\deqn{\widehat\Sigma_G=(MP_G-MP_{G:E})/(er),}
where \eqn{e} and \eqn{r} are the numbers of environments and replicates.
The same equations apply to the diagonal variances and off-diagonal
covariances. All raw component matrices are retained, including negative
variances or indefinite matrices. \code{object$moments$positive.semidefinite}
reports their numerical semidefiniteness; no matrix is repaired or truncated.

\code{st.lmm.components(object, "sscp")$strata} and
\code{object$moments$mean.products} are named lists of trait matrices for
every source and \code{Residual}.
\code{object$moments$strata} gives source DF;
\code{object$moments$EMS} contains expected mean-product coefficients;
\code{object$moments$moment.coefficients} contains the square system used
to estimate the components. The estimates are available as
\code{object$engine$sigma}, where the ordinary residual component is
named \code{units}. The moments-specific diagnostic record
\code{object$moments$covariance.components} retains the stratum name
\code{Residual}. Printing shows the strata, mean products,
coefficient matrix and all covariance component matrices automatically.
\code{anova} returns cached trait-wise fixed or requested random F tests.
\code{object$moments$stratum.geometry} stores the coefficient matrices
that weight each mean-product stratum in the fixed-coefficient covariance.
\code{st.lsm} uses this information for within-trait confidence intervals
and tests. Joint tests in \code{st.lsm} concern one trait per hypothesis.
Multivariate prediction and random-effect BLUPs use likelihood fits.
}
For a marker dataset, \code{y} is a snapshot of \code{gs.get.y(data.set, names = TRUE)} and \code{gen} is a factor of its individual identifiers, in marker-data order. Missing performance values are retained as \code{NA} and omitted from estimation by the ordinary model machinery. With \code{data="maize", ucov=list(gen="A")}, the quoted \code{"A"} requests the internally constructed VanRaden relationship matrix, with \code{coding="0 1 2"}. It is constructed on demand and reused within the call. All marker individuals enter the relationship matrix, including those without observed responses. No caller object is created or changed.

Only quoted strings are symbolic covariance requests. Unquoted \code{A} evaluates the actual R object named \code{A}, including when \code{data} names a stored data set. Prepared lists such as \code{U <- list(gen="A")} work identically to inline lists. Symbolic covariance requests require a stored marker or hybrid data-set input; ordinary data frames require actual numeric covariance matrices.

For stored hybrid data, \code{ucov=list(mother="A", father="A", "mother:father"="kro(a,a)")} requests parental GCA matrices and their mapped SCA product. The list names correspond to the actual parental factor names. The quoted product \code{"kro(A,A)"} is equivalent; the placeholders refer to the two parents in term order. The entire product is quoted, not only its operands. No R function named \code{kro} is needed. Other symbolic products are rejected; known general products require an actual \code{kronecker(M1,M2)} matrix. Strings are never evaluated as R code.

Successful converged top-level fits using quoted symbolic covariance requests also store fitted prediction state in their source marker datasets. Hybrid parental sources share that state. It is accessible through \code{\link{st.marker.results}} and \code{\link{st.predict}}. Failed fits and internal refits do not replace it. Numeric covariance matrices do not establish a new marker association. Simulation populations do not receive this state.

This interface fits the genomic second stage when stored performance values are genotype means or adjusted means. It does not fit a field-trial first stage or calculate first-stage weights. \code{st.load.performance.data} loads those responses, and \code{gs.set.performance.data} replaces them. Without \code{ucov}, genotype effects remain independent and no relationship matrix is constructed. The fitted object's \code{data} and \code{ucov} retain the actual model data and covariance matrices, so subsequent dataset changes do not alter an existing fit.

Fixed formulas support \code{+}, \code{-}, \code{:}, \code{*}, \code{/}, powers, \code{0}/\code{1}, parentheses and numeric expressions such as \code{I(x^2)}. The formula compiler walks R's parsed formula objects; it does not parse deparsed strings. Numeric transformations are evaluated in the formula/data environment. In a fixed formula, \code{.} expands to the other named data columns after
excluding response columns. In mixed models, grouping-only identifiers and
residual-structure columns are excluded from this automatic expansion;
explicitly naming such a column in the fixed formula still includes it.
Fixed-formula \code{offset(...)} terms and the optional argument offset
contribute known values to the response mean. They are summed once, retained
through refitting and updates, and reused for prediction according to the
offset rules above. Unordered factor labels that are all finite numbers are
sorted numerically; ordered levels retain their declared order. Factor
coding follows the explicit \code{contrasts} argument, then a factor's
contrasts attribute, then the global contrasts option, including its
separate setting for ordered factors. The fitted coding is retained for
prediction. Type III hypotheses use sum coding internally. Type II
containment requires matching continuous variables and nesting of
categorical variables.

For a multivariate response, genetic trait covariances for a random
\code{gen} term are available as \code{m$engine$sigma$gen}.
\code{cov2cor(as.matrix(m$engine$sigma$gen))} returns the corresponding
genetic correlation matrix. These are correlations between random genetic
effects; they are distinct from correlations between observed traits.

For a multivariate response, every bare random term and \code{units} defaults to an unstructured trait covariance, equivalent to \code{us(trait)}. \code{diag(trait):gen} and \code{diag(trait):units} specify independent traits. Each \code{at(env):gen} stratum receives its own block; \code{at(env,"local"):gen} selects levels. A multivariate random regression has a full covariance over trait by regression coefficient unless restricted explicitly.

For \code{random = ~gen} with \code{ucov = list(gen = A)}, the genetic effects have covariance \eqn{\sigma_g^2 A}; the variance multiplier is estimated by ML or REML. The supplied matrix is the relationship matrix itself; an inverse is not expected. No centring or rescaling is applied to the supplied matrix, so the reported variance component is a multiplier of its original scale. Singular genomic relationship matrices are supported. For multivariate effects or random regressions, the covariance between the coefficient vectors of groups \eqn{g} and \eqn{h} is \eqn{A_{gh}G}, where \eqn{G} is the estimated trait/regression covariance.

Names in \code{ucov} match the grouping part of a random term: \code{gen} also applies to \code{~Time | gen}, \code{us(trait):gen} and each independent stratum of \code{at(env):gen}. An interaction grouping term such as \code{gen:env} requires its own matrix and list name \code{"gen:env"}; a matrix named \code{gen} is not automatically expanded over environments. Residual covariance remains specified through \code{cov}.

All levels in a supplied matrix are retained in \code{object$engine$u}, including levels without phenotypic records. Conditional predictions for these levels use their relationship-based BLUPs. Levels absent from the supplied matrix are rejected because their relationships have not been specified; \code{allow.new.levels} does not supply missing relationships. The aligned matrices are retained in \code{object$ucov}; changing or removing the original matrix after fitting does not change predictions.

Missing response cells are retained as missing cells. A row with missing model predictors or grouping values is excluded for all its traits.
Without selection or observation weights, returned fitted values and residuals
are restored to the original row layout. Selection and weighted fits retain
the common analysis sample in the fitted object; the
\code{observation.rows} attribute maps those rows to the supplied data. Non-estimable fixed functions are returned as unavailable. Fitted values and residuals are matrices with one row per observation and one column per trait. Coefficient vectors are named.

The default requested fixed-effect inference is a Type III F-test with Satterthwaite denominator degrees of freedom. Certified balanced orthogonal terms use exact F ratios when available; other terms retain Satterthwaite inference, as in \pkg{lmerTest}. Satterthwaite inference uses the observed covariance-parameter information and derivatives of fixed-effect covariance. Tests are computed separately by trait and jointly across traits. Zero independent random-variance components use the first-order standard-deviation information approximation; singular nonnegative information uses its positive-eigenvalue generalized inverse. These cases are identified in the engine information status. Negative-curvature or nonconvergence cases retain unavailable inference rather than substituting another degrees-of-freedom method. Asymptotic tests are Wald tests and display infinite denominator DF. Kenward--Roger uses expected REML information in linear covariance-component coordinates, adjusts the fixed-effect covariance, and applies hypothesis-specific denominator DF and F scaling. Nonconverged or unidentified covariance models do not yield Kenward--Roger inference. \code{st.lsm(object, L=L, joint=TRUE)} provides explicit hypotheses and the scaling factor. The raw engine covariance is retained for estimation; \code{vcov} and \code{st.lsm} expose the adjustment.

With \code{random.test="auto"}, a certified balanced orthogonal design uses exact random-effect F tests when every random term has one exact denominator; otherwise all random terms use LRT. Explicit \code{random.test="LRT"} removes each whole random term in turn and refits the reduced model with the same observations, fixed effects and estimation method. The statistic is twice the full-minus-reduced log-likelihood; the chi-squared degrees of freedom are the difference in the number of free covariance parameters. REML fits use REML likelihoods. This corresponds to whole-term deletion in \code{lmerTest::ranova(reduce.terms=FALSE)}, including joint deletion of all strata and coefficients of a random regression or multivariate term. It does not perform stepwise slope reductions. The chi-squared reference distribution is approximate, without a boundary-mixture correction; the optional bootstrap provides an alternative when that approximation is unsuitable. LRT tests use no simulations and do not consume random numbers. A failed reduced fit or an invalid likelihood comparison gives unavailable inference with a reason.

Bootstrap random-term tests refit full and reduced models with identical fixed effects, using the fitted method (ML or REML), simulated under the reduced model. They use R's random-number generator; \code{set.seed} controls repeatability. All strata and trait/regression coefficients belonging to a random term are tested together. The p-value is \code{(exceedances+1)/(replicates+1)}. Monte Carlo standard error and refit failures are reported; any failed refit causes the p-value to be withheld. Bootstrap may be expensive: \code{random.test="none"} is useful for exploratory fitting. Printing never simulates or refits.

Marginal and conditional R-squared use the sample variance of fitted fixed effects and the observation-averaged diagonal of random/residual covariance, separately by trait over observed responses. RMSE uses conditional residuals. \code{AIC} and \code{BIC} dispatch through the log-likelihood method; the counted parameters are the estimable fixed coefficients plus estimated covariance parameters. Usual restrictions on comparing REML likelihoods across different fixed models apply.

New fits retain fixed-effect and random-regression prediction recipes for
orthogonal polynomials, scaling and spline bases, including their training
coefficients, center, scale and knots. Prediction and least-squares means reuse
these recipes. Transformations are calibrated before row subsetting and
missing-value removal. Older saved fits without the required metadata must
be refitted or use precomputed columns. Formulas may call user functions,
so their environments must remain available to evaluate new data. All
fitted values, tests and statistics can be serialized without external
pointers.

\subsection{Prediction of selected random effects and uncertainty}{
The following extensions apply to ML and REML fits. The default
\code{predict(object)$fit} contains the marginal fitted means, including
offsets. Conditional
prediction adds the selected fitted random coefficients. Random-only
prediction returns those random contributions and excludes fixed effects
and offsets. Formula selectors can select an intercept or slope within a
correlated block. Selection changes the prediction target only: all fitted
random components remain in the observation covariance used for conditioning.
A named-block selector therefore differs from deleting terms and refitting.

The target fixed design is \eqn{X_*}; \eqn{C_*} is the prior covariance of
its selected random contribution, \eqn{K_*} is its covariance with the observed
responses, and \eqn{V} is the full observed-response covariance. With estimated
fixed effects, the latent-mean prediction-error covariance is
\deqn{C_* - K_* V^{-1} K_*' + D\,\mathrm{Var}(\widehat\beta)D',
       \qquad D = X_* - K_* V^{-1} X.}
Random-only targets use \eqn{X_*=0}; marginal targets use
\eqn{C_*=K_*=0}. \code{se.fit} returns the square roots of the diagonal.
This includes the joint fixed/random prediction-error correction and is
not obtained by adding separate coefficient standard errors. Covariance
parameters are held fixed, future observation residuals are excluded, and
no covariance-parameter uncertainty is implied. Nonestimable fixed
combinations are returned as unavailable.

For independent unseen groups allowed by \code{allow.new.levels=TRUE}, the
random mean is zero and its prior variance contributes to prediction error.
For named levels already present in a supplied relationship matrix,
relationship-based prediction and uncertainty use the covariance with the
observed levels, even when the level had no observed response. A new
covariance stratum has no estimated variance, so its conditional standard
error is unavailable and requesting one produces an error. New data for a
direct observation kernel require an explicitly defined target kernel and
are not supported by this interface.

\code{st.lmm.components(object, c("fixed", "covariance"))} provides the
physical parameter layout accepted by \code{newparams}. Supplying covariance
parameters performs one covariance evaluation and GLS/BLUP solve without
optimization. If \code{fixed} is omitted, fixed coefficients are estimated
at the supplied covariance and their estimation contributes to prediction
error. If \code{fixed} is supplied, its values are treated as known and the
random conditional mean is recomputed from the residuals relative to those
values. Its conditional covariance is \eqn{C_* - K_*V^{-1}K_*'}; the
fixed-effect estimation term is omitted. Accordingly, supplying the fitted
coefficient values explicitly gives the same point prediction but can give
a different standard error. Supplying complete fixed coefficients also
defines predictions involving columns aliased during fitting. The fitted
object is not modified.

The returned uncertainty list has attributes describing its interpretation,
whether fixed coefficients are treated as known, the selected coefficients,
and the exclusion of future residual noise. With \code{na.exclude},
training-row predictions are padded at omitted positions; missing-response
prediction rows retained under the default policy receive predictions.
The \code{re.form} and \code{newparams} extensions are not available for
\code{method="moments"}; its existing marginal/conditional prediction and
interval interface is unchanged.
}
}
\section{Covariance formula syntax}{
The following formulas are interpreted by the \code{st.lmm} model
compiler. The symbols \code{trait} and \code{units} are reserved in these
covariance terms; no columns with these names need to be added to \code{data}.
\code{units} identifies the original data rows. In terms containing
\code{units}, different rows have independent residuals, while traits
measured on the same row may have residual covariances. The separate
\code{occasion | group} form described below permits residual correlation
between repeated observations in different rows.

\subsection{Trait covariance and independent strata}{
\tabular{ll}{
\strong{Term} \tab \strong{Meaning} \cr
\code{us(trait)} \tab Unstructured covariance across traits. \cr
\code{diag(trait)} \tab A separate variance per trait; zero trait covariances. \cr
\code{homdiag(trait)} \tab A common residual variance; zero trait covariances. \cr
\code{cs(trait, hom=TRUE)} \tab A common residual variance and correlation. \cr
\code{cs(trait, hom=FALSE)} \tab Separate variances and a common correlation. \cr
\code{ar1(trait)} \tab A common residual variance and AR(1) correlations. \cr
\code{ar1(trait, hom=FALSE)} \tab Separate variances and AR(1) correlations. \cr
\code{at(loc)} \tab A separate, independent block for each location. \cr
\code{at(loc, "HOH")} \tab Only the block for location HOH. \cr
\code{at(loc, c("HOH", "EWE"))} \tab Two separate location blocks. \cr
}

For a matrix response such as \code{cbind(GYLD, TKW, PLH)}, omitting the
trait modifier means \code{us(trait)}, both for random effects and residuals.
Thus \code{~units} and \code{~us(trait):units} specify the same residual
model. An unstructured three-trait block has three variances and three
covariances; it does not constrain the traits to have the same variance.
\code{diag(trait)} also gives each trait its own variance, but fixes all
off-diagonal covariances to zero. These restrictions apply to the component
in that term, not to other random or residual components.

\code{at(loc)} changes which parameters are shared between locations;
\code{diag(trait)} changes which traits may covary. Combining them as
\code{diag(trait):at(loc):units} estimates one residual variance for every
trait by location combination. Selecting several levels with \code{at()}
still gives a separate block for each selected level; it does not pool them
into a common variance group.
}

\subsection{Single-trait residual formulas}{
For a single numeric response, each of the following independent-residual
blocks has one variance.
The counts \eqn{L} and \eqn{E} are the numbers of observed locations and
environments, respectively, with no missing grouping values.
\tabular{lll}{
\strong{Formula for cov} \tab \strong{Residual variances} \tab \strong{Count} \cr
\code{~units} \tab One common variance. \tab 1 \cr
\code{~at(loc):units} \tab One per location. \tab \eqn{L} \cr
\code{~at(env):units} \tab One per environment. \tab \eqn{E} \cr
}
The meaning of \code{env} depends on how that variable was constructed;
its name alone does not define an environment or imply a location.
}

\subsection{Correlated residuals for repeated observations}{
For a univariate ML or REML fit, \code{cov = ~ar1(0 + Visit) | Subject}
defines AR(1) residual covariance between visits within each subject.
Different subjects have independent residuals. \code{Visit} must be an R
factor; its declared levels, including unobserved levels, determine the
occasion order. Both ordered and ordinary factors retain that order.
The coefficient formula must consist of this single factor with no
intercept, as in \code{0 + Visit} or \code{Visit - 1}. Numeric regressions,
character vectors, additional predictors and contrast bases are rejected.
Each residual stratum may have at most one observed row per subject and
visit. Missing responses and incomplete sets of visits are allowed.

With \eqn{q} occasion levels, the following modifiers are supported on the
left of \code{| Subject}:
\tabular{lll}{
\strong{Modifier} \tab \strong{Residual covariance} \tab \strong{Parameters} \cr
\code{us(0 + Visit)} \tab Unstructured. \tab \eqn{q(q+1)/2} \cr
\code{diag(0 + Visit)} \tab Separate variances; zero covariances. \tab \eqn{q} \cr
\code{homdiag(0 + Visit)} \tab Common variance; zero covariances. \tab 1 \cr
\code{cs(0 + Visit, hom=TRUE)} \tab Common variance and correlation. \tab 2 \cr
\code{cs(0 + Visit, hom=FALSE)} \tab Variances differ; correlation shared. \tab \eqn{q+1} \cr
\code{ar1(0 + Visit)} \tab Common variance; AR(1) correlations. \tab 2 \cr
\code{ar1(0 + Visit, hom=FALSE)} \tab Separate variances; AR(1) correlations. \tab \eqn{q+1} \cr
}
With one occasion, each structure has one variance and no correlation
parameter. The default is \code{hom=TRUE} for \code{ar1()} and
\code{hom=FALSE} for \code{cs()}. For visits at positions \eqn{i} and
\eqn{j}, AR(1) specifies residual covariance
\eqn{s_i s_j\rho^{|i-j|}}; compound symmetry replaces the correlation by
\eqn{\rho} for distinct positions. Homogeneous structures use a common
\eqn{s_i}. Correlation bounds and variance constraints are the same as in
\sQuote{Random-coefficient covariance structures}.

AR(1) lags refer to factor-level positions, not differences between numeric
labels. Intermediate occasions in the factor levels represent intended gaps. For example, levels \code{c("1", "2", "3")} retain a lag of two
between occasions 1 and 3 even if occasion 2 is unobserved. Homogeneous
structures can retain entirely unobserved occasion levels. A separately
estimated occasion variance requires observed information in that stratum;
otherwise the fit is rejected. A missing level is never silently dropped
to shorten AR(1) distances. Insufficient information for a correlation is
reported by the covariance-information diagnostics. When all observed
within-group pairs have even lags, the likelihood cannot distinguish
positive from negative AR(1) correlation; its sign is unidentified, and
the profile inventory does not offer an ordinary interval for signed
\eqn{\rho}.
If some heterogeneous AR(1) variances are zero and no observed pair has two
positive standard deviations, the correlation can be unidentified while
still affecting directions that increase a zero variance. The optimizer
checks these directions. If the remaining boundary curvature cannot be
certified, it reports nonconvergence rather than accepting that boundary.

The grouping side can include interactions and independent strata, for
example \code{~ar1(0 + Visit) | at(Treatment):Subject}, which estimates a
separate covariance per treatment. Lists add residual components. The
existing \code{~at(Treatment):units} notation continues to describe
independent residuals with a treatment-specific variance.
With precision weights, covariance between rows \eqn{i} and \eqn{j} is
divided by \eqn{\sqrt{w_iw_j}}. Random effects and \code{ucov} may be used
alongside these residual structures.

Refitting, simulation and construction of future residual covariance retain
the fitted occasion levels and their order; an unknown occasion level is
an error. Mean predictions retain their documented fixed/random-effect
meaning; they do not add a fitted residual process to the predicted mean.
These structures require ML or REML. The inference restrictions described
for the corresponding covariance families also apply to residual blocks.
Correlated residuals use the exact block or dense covariance backend;
forcing an incompatible sparse backend is an error. Model comparisons
report likelihood and information criteria, but do not infer an LRT
reference distribution for changes in the residual covariance family.
}

\subsection{Separable spatial residual covariance}{
For a univariate ML or REML fit,
\code{cov = ~ar1(Row):ar1(Column)} estimates a homogeneous separable
AR(1) by AR(1) plot-error covariance. The colon denotes a direct product,
using the formula style of ASReml-R 3. The existing \code{cov} argument
continues to specify the residual covariance.
There are three estimated parameters: one residual variance and one
correlation for each factor. For two plots in rows \eqn{r,r'} and columns
\eqn{c,c'}, their residual covariance is
\deqn{\sigma_e^2\rho_r^{|r-r'|}\rho_c^{|c-c'|}.}
Both variables must be distinct R factors, each with at least two declared
levels. Declared factor levels determine the positions and spacing, including
unobserved intermediate rows or columns. Data rows may be in any order;
the first factor varies fastest in the retained coefficient covariance
matrix. For the formula above, rows within columns give
\deqn{R=\sigma_e^2\{A_C(\rho_c)\otimes A_R(\rho_r)\},}
where \eqn{[A_m(\rho)]_{ij}=\rho^{|i-j|}}.
Missing responses and incomplete grids are supported. A grid may have at
most one observed response per row-column combination. The form
\code{~ar1(Row):ar1(Column) | Field} gives independent fields sharing the
same spatial parameters; \code{at()} strata on the grouping side estimate
separate covariance blocks. Factor levels retain their order for refitting
and future residual covariance construction. Unknown coordinate levels
are rejected.

Random replicate and incomplete-block terms may be retained, for example
\code{random = ~Rep + RowBlk + ColBlk}. Precision weights divide residual
covariances by \eqn{\sqrt{w_iw_j}}. A list of residual formulas adds
components; it does not multiply them. In particular,
\code{list(~ar1(0 + Row) | Column, ~ar1(0 + Column) | Row)}
is an additive row-column covariance, not a separable AR(1) by AR(1)
covariance.

Analytic covariance derivatives are used for likelihood optimization and
Satterthwaite inference. Asymptotic inference is also available.
Kenward--Roger inference and likelihood profiles are not implemented for
the separable product. Correlations have bounds \eqn{[-1,1]}; singular
parameter combinations are admissible only if the full model covariance
remains positive definite on the observed responses. With observations
only at even lags in a direction, the corresponding correlation sign is
unidentified. The optimizer checks nonzero-correlation alternatives before
accepting a zero-correlation stationary fit. A zero spatial variance whose
two-dimensional covariance cone cannot be certified is reported as
nonconverged, rather than accepted as a verified boundary optimum.
The printed summary shows the variance and both correlations; the full
coefficient covariance remains available from \code{st.lmm.components()}.
}

\subsection{Multivariate residual formulas}{
For \eqn{T} response traits, \eqn{U=T(T+1)/2} denotes the number of
parameters in one unstructured trait covariance matrix. With all traits
observed in every stratum, the following are the numbers of free residual
covariance parameters; random-effect parameters are additional.
\tabular{lll}{
\strong{Formula for cov} \tab \strong{Residual structure} \tab \strong{Count} \cr
\code{~units} \tab Unstructured trait covariance. \tab \eqn{U} \cr
\code{~us(trait):units} \tab Same unstructured model. \tab \eqn{U} \cr
\code{~diag(trait):units} \tab Separate trait variances. \tab \eqn{T} \cr
\code{~homdiag(trait):units} \tab Common variance; zero covariance. \tab 1 \cr
\code{~cs(trait, hom=TRUE):units} \tab Common variance and correlation. \tab 2 \cr
\code{~cs(trait, hom=FALSE):units} \tab Variances differ; correlation shared. \tab \eqn{T+1} \cr
\code{~ar1(trait):units} \tab Common variance; AR(1). \tab 2 \cr
\code{~ar1(trait, hom=FALSE):units} \tab Separate variances; AR(1). \tab \eqn{T+1} \cr
\code{~at(loc):units} \tab Unstructured, per location. \tab \eqn{LU} \cr
\code{~us(trait):at(loc):units} \tab Same location model. \tab \eqn{LU} \cr
\code{~diag(trait):at(loc):units} \tab Trait variances per location. \tab \eqn{LT} \cr
\code{~at(env):units} \tab Unstructured, per environment. \tab \eqn{EU} \cr
\code{~diag(trait):at(env):units} \tab Trait variances per environment. \tab \eqn{ET} \cr
}

For example, \code{~diag(trait):at(loc):units} assumes, for unit \eqn{i}
at location \eqn{l},
\deqn{\mathrm{Var}(e_i) =
  \mathrm{diag}(\sigma^2_{1l},\ldots,\sigma^2_{Tl}).}{
  Var(e_i) = diag(sigma^2_(1,l), ..., sigma^2_(T,l)).}
Using \code{~at(loc):units} instead replaces that diagonal matrix by a
location-specific unstructured matrix. In both cases, different data rows
have zero residual covariance. Random effects can still induce covariance
between their observations.

The \code{homdiag}, \code{cs} and \code{ar1} trait modifiers require ML
or REML and may also be combined with \code{at()} strata. Their parameter
counts above assume at least two traits; with one trait there is one
variance and no correlation parameter. For AR(1), the left-to-right
response order in \code{cbind(y1, y2, y3)} defines positions 1, 2 and 3.
Changing that order can change the covariance model. Traits on different
data rows remain residual-independent. These formulas do not combine a
trait covariance with an additional repeated-occasion covariance.
}

\subsection{Location versus year by location: the barley MET example}{
The \code{barley.met} dataset contains five recorded traits, five
locations and three years. The following covariance examples use the three
responses \code{GYLD}, \code{TKW} and \code{PLH}. Here \code{env}
is a year-by-location factor:
\preformatted{
data("v-lmm-barley-met", package = "SelectionTools")
barley.met$env <- interaction(barley.met$year, barley.met$loc,
                              drop = TRUE)
}
There are 15 observed environments. Therefore:
\tabular{lll}{
\strong{Formula for cov} \tab \strong{Residual structure} \tab \strong{Count} \cr
\code{~at(env):units} \tab Full trait covariance per year by location. \tab 90 \cr
\code{~diag(trait):at(env):units} \tab Trait variances per year by location. \tab 45 \cr
\code{~at(loc):units} \tab Full trait covariance per location. \tab 30 \cr
\code{~diag(trait):at(loc):units} \tab Trait variances per location. \tab 15 \cr
}

The location formulas share each location's residual parameters across
years. The environment formulas allow them to differ between years within
the same location. Equivalently, for these data,
\code{~at(year):at(loc):units} specifies separate unstructured blocks for
the observed year by location combinations, and
\code{~diag(trait):at(year):at(loc):units} specifies separate variances
for each trait in each such combination.

With eight bare random terms, each term has an unstructured three-trait
covariance matrix, giving 48 random-effect parameters in total.
Together with \code{~at(env):units}, this gives 138 covariance
parameters. Using \code{~diag(trait):at(loc):units} gives 63 instead.
This is a different residual model; reducing the parameter count does not
by itself establish convergence or justify the covariance restriction.
}

\subsection{Selected levels and combinations of formulas}{
An \code{at()} level selector activates a covariance component only for
the named levels; it does not subset the data. For example, residual terms
using only \code{at(loc, "HOH")} do not provide residual variation for
observations at the other locations. A fit for HOH alone requires a
corresponding data subset; a fit including other locations requires residual
terms covering them.
The following level labels provide a complete partition of the barley MET
locations:
\preformatted{
cov = list(
    ~diag(trait):at(loc, c("HOH", "EWE")):units,
    ~us(trait):at(loc, c("EIN", "OLI", "HFH")):units
)
}
This estimates three independent trait variances at each of HOH and
EWE, and a full three-trait covariance matrix separately at each of the
remaining locations: 24 residual parameters in total. A list adds the
specified covariance components. Overlapping residual terms require
an intended additive covariance model with identifiable parameters.
}

\subsection{The same modifiers for random effects}{
For a multivariate response, the corresponding random formulas are:
\tabular{ll}{
\strong{Formula for random} \tab \strong{Random-effect structure} \cr
\code{~gen} \tab One unstructured trait matrix for genotype effects. \cr
\code{~us(trait):gen} \tab The same genotype model, stated explicitly. \cr
\code{~diag(trait):gen} \tab Separate genotype variances per trait; zero trait covariances. \cr
\code{~at(loc):gen} \tab Independent genotype effects per location; a trait matrix per location. \cr
\code{~diag(trait):at(loc):gen} \tab Independent trait by location genotype components. \cr
\code{~gen + gen:loc} \tab Common genotype effects plus independent genotype by location effects. \cr
}
The terms \code{gen:loc} and \code{at(loc):gen} differ: the former shares
one trait covariance matrix across genotype by location effects, whereas
the latter estimates a separate matrix per location. Adding \code{gen}
also induces genetic covariance between observations of the same genotype
at different locations. Without a matching \code{ucov} entry, distinct
grouping levels are independent; with \code{ucov}, the supplied grouping
covariance is combined with the trait covariance as described in Details.
}

The residual \code{at(...):units} notation resembles ASReml-R 3 syntax.
SelectionTools also supports the direct-sum subset documented below;
this does not implement every ASReml-R 4 \code{dsum} option.
}
\section{Covariance composition}{
For univariate ML/REML models, colon-separated factor covariance modifiers
specify separable products. Each axis must be a distinct, directly named R
factor with at least two declared levels. Supported product factors are
\code{ar1(Factor)}, \code{ar1(Factor, hom=TRUE)},
\code{cs(Factor, hom=TRUE)} and \code{homdiag(Factor)}.
A product has one variance multiplier, not one variance per factor.
Each AR1 or CS factor adds one correlation; homdiag adds none.
CS correlation bounds depend on that axis's number of levels.

\tabular{ll}{
\strong{Request} \tab \strong{Meaning} \cr
\code{cov=~ar1(Row):cs(Column,hom=TRUE) | Field}
\tab Independent fields with shared product parameters. \cr
\code{cov=list(~ar1(Row):cs(Column,hom=TRUE) | Field, ~units)}
\tab Product covariance plus an independently estimated nugget. \cr
\code{cov=~(ar1(Row):cs(Column,hom=TRUE) | Field)+units}
\tab Equivalent additive formula; parentheses delimit the bar term. \cr
\code{cov=~dsum(~ar1(Row):cs(Column,hom=TRUE) | Env, group=Field)}
\tab Environment-specific product parameters, shared across fields within
an environment; all field blocks are independent. \cr
\code{random=~ar1(Visit):homdiag(Season) | Subject}
\tab A product covariance among random coefficients within subjects. \cr
}

Products can have more than two axes and can be parenthesized associatively.
The first axis varies fastest in the retained coefficient covariance matrix.
AR1 lags use positions in declared factor levels, including unobserved levels.
They are not numeric distances between observation labels. Missing grid cells
are supported; their declared levels determine the spacing.
Residual products require at most one observed response per cell and group;
random coefficients can have repeated measurements in a cell.
A named \code{ucov} matrix combines the coefficient product with the known
covariance between grouping levels. Precision weights scale the whole
observation covariance by the reciprocal square roots of the weights.

\code{dsum} accepts exactly one unnamed one-sided formula
\code{dsum(~covariance | stratum)} and optional \code{group=Subject}.
The stratum and optional group must be distinct directly named data factors.
For example, \code{dsum(~units | Env)} equals \code{at(Env):units}, and
\code{dsum(~ar1(0+Visit) | Treatment, group=Subject)} equals
\code{ar1(0+Visit) | at(Treatment):Subject}. Existing within-row multivariate
trait covariance also works with \code{dsum(~units | Env)}.
The constructor is residual-only. Independent random strata are specified by
\code{at(stratum)} in the grouping expression.
Unknown named arguments, nested dsum, repeated axes and malformed formulas
are errors. A list of formulas specifies additive covariance components.

Heterogeneous product factors (\code{hom=FALSE}, \code{diag}) and arbitrary
unstructured product factors are not implemented. The supported
\code{us(trait)} exception for multivariate repeated records is described
below. CS product factors must explicitly request \code{hom=TRUE}; their
ordinary non-product default is unchanged. The fitted covariance layout,
axis names, bounds and parameters are retained for refitting, prediction,
simulation and \code{st.lmm.objective}. User-supplied physical covariance
matrices must obey the fitted product shape, not merely be positive
semidefinite. Separable families support Satterthwaite and asymptotic
inference; Kenward--Roger and likelihood profiles are rejected explicitly.
At a zero product variance with unidentified correlations, convergence is
not certified unless its covariance cone can be checked. Existing
\code{ar1(Row):ar1(Column)} residual models retain their established path.
}
\section{Multivariate repeated records}{
For a matrix response with at least two traits, one unstructured trait
covariance can be combined with one occasion-correlation factor inside an
independent residual group:
\preformatted{
st.lmm(cbind(y1, y2) ~ treatment + x,
    cov = ~us(trait):ar1(Visit) | Subject,
    data = dat, ddf = "asymptotic", random.test = "none")

cov = ~us(trait):cs(Visit, hom=TRUE) | Subject
cov = ~us(trait):homdiag(Visit) | Subject
}
Here \code{trait} is the virtual response-trait dimension; it is not a
required data column. \code{Visit} must be an R factor with at least two
declared levels. \code{Subject} is one directly named grouping factor.
Its distinct levels define independent residual groups that share all
trait and occasion parameters. Reversing the two product factors gives
the same canonical model and coefficient order.

With \eqn{p} traits and \eqn{q} declared occasions, the coefficient
covariance is \eqn{\Sigma_T\otimes R_O}. Occasion varies fastest within
trait in the returned coefficient covariance. \eqn{\Sigma_T} is one full
\eqn{p\times p} trait covariance matrix. \eqn{R_O} has diagonal one and
is AR1, homogeneous compound symmetry, or identity for \code{homdiag}.
The number of estimated covariance parameters is \eqn{p(p+1)/2+1} for
AR1 or CS and \eqn{p(p+1)/2} for homdiag. There is no second occasion
variance multiplier. In particular, homdiag occasions are independent,
while traits measured on the same occasion can remain correlated.
For AR1, correlation between occasions at declared level positions
\eqn{a,b} is \eqn{\rho^{|a-b|}}. For CS, all distinct occasions have
correlation \eqn{\rho}, whose lower bound is \eqn{-1/(q-1)}.

Missing trait responses can be retained independently under the default
missing-response handling or \code{na.action=na.pass}. A missing cell does
not require dropping the other observed traits in its row. Whole missing
occasions and incomplete subject sequences are supported. The declared
factor levels determine the occasion spacing: observed levels 1 and 3
have AR1 correlation \eqn{\rho^2}
even when level 2 has no data row. Explicit row-removing actions such as
\code{na.omit} retain their ordinary whole-row semantics. There may be at
most one observed response per subject, occasion and trait; complementary
missing-trait rows at the same subject/occasion key are allowed.

The observed cells must identify the requested parameters. The compiled
check tests covariance directions at generic interior parameter values,
using the fixed-effect projection for REML. Missing cross-trait support,
fixed-effect absorption and dependence between random and residual
components can cause an explicit error. For example, an unrestricted
subject-level trait random intercept together with a CS trait-product
residual generally cannot identify the additional occasion correlation.
This check concerns local covariance identifiability, rather than a proof
of unique global parameter labels in every possible nonlinear mixture.

To certify an AR1 correlation sign, this implementation additionally
requires at least one observed odd-lag pair within a trait and subject.
An even-lag-only sequence cannot distinguish positive and negative rho.
Patterns that would determine the sign only through cross-trait covariance
are currently rejected as uncertified. Odd gaps larger than one are
supported; fitting probes nonzero correlations to escape a flat starting
point at zero. At a locally flat or otherwise rank-deficient information
matrix, covariance-parameter uncertainty must be treated as unavailable.

Fitting uses ML or REML and compiled L-BFGS-B. Explicit AI requests are
rejected. The default optimizer comparison in \code{st.lmm.allfit} uses
L-BFGS-B for this family. Warm/cold refits, response-unit changes, update,
serialization, ordinary prediction and simulation retain the trait/occasion
layout. The fitted covariance is exported on the original response scale,
including when occasion frequencies differ or a declared occasion is absent.
\code{st.lmm.objective} reports trait log-SDs and transformed partial
correlations, followed by \code{rho.occasion} where applicable. Its trait
coordinate row/column indices refer to the trait dimensions; SD scale
multipliers convert back to original response units. These optimizer
coordinates differ from a vector of physical marginal trait correlations.

Supplied covariance matrices for prediction or simulation must have the
full trait-by-occasion dimensions, matching coefficient names and the
requested product form. Positive semidefiniteness alone is insufficient.
A positive-semidefinite singular product can be simulated; a fitted singular
or numerically unresolved boundary does not certify ordinary likelihood
inference. Satterthwaite inference requires usable covariance-parameter
information; asymptotic fixed-effect inference conditions on the fitted
covariance. Kenward--Roger adjustment and covariance profiles are currently
unsupported and are rejected explicitly.

This extension is residual-only and contains one virtual \code{us(trait)}
and one ordinary occasion axis per product. Additional axes, random
trait-products, other trait modifiers within products, heterogeneous
occasion scales, nested/stratified grouping expressions and direct sums
of trait-products are not yet supported. Ordinary supported random-effect
terms can accompany an identifiable trait-product residual. Explicit
fixed/shared covariance constraints remain unavailable for models containing
this estimated US trait family. Multivariate precision weights and offsets
retain their existing unsupported-interface errors. Multivariate inheritance
from earlier fits is a separate extension and is not implied by this syntax.
}
\section{Fixed and shared covariance parameters}{
\code{constraints} specifies parameter values in original response and
coefficient units, independently of internal response or design scaling.
Constraints require \code{method="ML"} or \code{method="REML"};
\code{method="moments"} is not supported.
The named list accepts exactly \code{fixed} and \code{shared}. Selectors
must come from the compiled model's inventory:
\preformatted{
inventory <- st.lmm.objective(fit)$constraints
inventory$parameter
inventory$variance
}
Canonical selectors begin with \code{sd:random:}, \code{sd:residual:},
\code{cor:random:} or \code{cor:residual:}, followed by the exact block and
coefficient names. A homogeneous variance uses \code{common}; a structured
correlation uses \code{rho}; separable products use \code{rho.1},
\code{rho.2}, and the corresponding axis numbers. The \code{variance}
inventory provides \code{var:} aliases for variance values. For example,
\code{var:random:gen:(Intercept)} fixes a variance, whereas
\code{sd:random:gen:(Intercept)} fixes its standard deviation. Inventory
\code{value} entries use the canonical SD or correlation units, not the
variance-alias units. Block names include any retained stratum labels;
the fitted names identify the valid selectors.

\preformatted{
st.lmm(y ~ 1, random = ~gen, data = dat,
    constraints = list(fixed = c("var:random:gen:(Intercept)" = 2)),
    ddf = "asymptotic", random.test = "none")

st.lmm(y ~ 1, random = ~gen + env, data = dat,
    constraints = list(shared = list(c(
        "sd:random:gen:(Intercept)", "sd:random:env:(Intercept)"))),
    ddf = "asymptotic", random.test = "none")
}

Every \code{shared} entry contains at least two selectors of the same
physical meaning. Each group consistently uses SD selectors, variance selectors
or correlation selectors. A parameter
cannot belong to more than one group. A fixed member fixes its whole shared
group; multiple fixed values in that group must agree in physical units.
Unknown or ambiguous names, duplicated declarations, nonfinite values,
incompatible units and infeasible parameter domains are errors. Correlation
groups use the intersection of their compiled bounds, including dimension-
dependent compound-symmetry bounds. Random variances may be fixed at zero
where the full marginal covariance remains valid.

Supported estimated blocks are diagonal, homogeneous diagonal, compound
symmetry, AR1 and the supported separable products. Models containing an
estimated unstructured covariance block, including a US trait-by-occasion
product, are currently rejected when explicit constraints are requested. The ordinary unconstrained unstructured interface
is unchanged. Fixed known covariance matrices remain fixed inputs.

Constrained fitting optimizes only independent parameters using compiled
L-BFGS-B. \code{optimizer="auto"} selects that path; an explicit AI request
is rejected. Starts are reconciled with the fixed/shared declarations before
optimization. Explicit objective parameter evaluations must already satisfy
the constraints and are never silently projected. With every covariance
parameter fixed, only the GLS/BLUP calculation is required.

The engine retains the full optimizer-coordinate vector for reconstruction
and reports its \code{full.count} and independent \code{free.count} in
\code{engine$constraints}. AIC, BIC and log-likelihood degrees of freedom count
independent estimated covariance parameters. The constraint inventory records
fixed/free/shared status and the independent coordinate index. Satterthwaite
inference uses covariance-parameter uncertainty projected through the reduced
parameter mapping; fixed parameters have zero parameter uncertainty.
The complete raw information matrix describes directions outside the
constraint surface as well and need not have the same definiteness or rank.
For constrained fits, \code{engine$uncertainty$free.information} and
\code{free.vcov} retain the reduced information and covariance matrices in
the \code{free.parameters} order when parameter uncertainty is computed.

Asymptotic and Satterthwaite fixed-effect inference, ordinary prediction and
simulation, serialization and response/weight refits retain constraints.
\code{random.test="none"} applies here. Kenward--Roger adjustment, covariance likelihood
profiles, automatic random-term deletion, subset selection and model-comparison
tests with explicit constraints are not yet implemented and fail explicitly.
Physical \code{newparams} overrides are rejected for constrained objects;
the intended covariance model requires an explicit fit. Optimizer comparisons use
L-BFGS-B and retain the same physical constraints for every supported start.
}
\section{Random-coefficient covariance structures}{
Both embedded \code{y ~ x + (1+x || group)} and separate
\code{random=~1+x || group} specify independent random coefficients,
equivalent to \code{random=~diag(1+x) | group}. Parenthesized additive
bar terms and lists of random formulas are supported consistently.
For factor predictors, independence applies to the expanded coefficient
columns in the chosen contrast basis, not to whole multi-column terms.
\code{||} cannot accompany covariance modifiers such as \code{cs()}
or \code{us()}; \code{|} supports these modifiers.
The operator \code{||} is not permitted in \code{cov}.

For a univariate ML or REML fit, modifiers around the coefficient formula
specify its covariance within each grouping level. For example,
\code{random = ~ar1(0 + year) | gen} uses one random coefficient for every
level of \code{year}, with an AR1 covariance among those coefficients for
each genotype. These modifiers do not change the fixed-effects formula.

\tabular{lll}{
\strong{Modifier} \tab \strong{Covariance} \tab \strong{Parameters} \cr
\code{us(...)} \tab Unstructured covariance. \tab \eqn{q(q+1)/2} \cr
\code{diag(...)} \tab Separate variances, zero covariances. \tab \eqn{q} \cr
\code{homdiag(...)} \tab One common variance, zero covariances. \tab 1 \cr
\code{cs(..., hom=TRUE)} \tab Common variance and common correlation. \tab 2 \cr
\code{cs(..., hom=FALSE)} \tab Separate variances and common correlation. \tab \eqn{q+1} \cr
\code{ar1(..., hom=TRUE)} \tab Common variance and AR1 correlations. \tab 2 \cr
\code{ar1(..., hom=FALSE)} \tab Separate variances and AR1 correlations. \tab \eqn{q+1} \cr
}
Here \eqn{q} is the number of random coefficient columns. With one column,
each structure estimates one variance and has no correlation parameter.
The default is \code{hom=FALSE} for \code{cs()} and \code{hom=TRUE} for
\code{ar1()}; \code{hom} must be a single literal \code{TRUE} or
\code{FALSE}. Unmodified random regressions and \code{us()} retain an
unstructured covariance. \code{diag()} continues to estimate a separate
variance for each coefficient.

For compound symmetry, \eqn{G_{ij}=s_i s_j \rho} for distinct coefficients
and \eqn{G_{ii}=s_i^2}, with
\eqn{-1/(q-1) \leq \rho \leq 1}. For AR1,
\eqn{G_{ij}=s_i s_j \rho^{|i-j|}}, with
\eqn{-1 \leq \rho \leq 1}. Homogeneous structures constrain all
\eqn{s_i} to the same nonnegative value. These constraints refer to the
original coefficient basis, before numerical column scaling.

AR1 distances are positions in the expanded coefficient basis. They are
not distances between numeric labels or observation times. For a factor,
\code{0 + year} specifies one coefficient per level. An ordered factor
states the intended order explicitly. An intercept, additional
numeric terms or different factor coding changes that basis and therefore
changes the covariance model. The coefficient names and their order are
retained on \code{object$engine$sigma}.

For random effects, these modifiers apply to coefficient blocks in
univariate likelihood fits. They can be combined with other random terms,
structured residual covariance and named matrices in \code{ucov}.
Their residual uses are described in \sQuote{Covariance formula syntax}.
They are not supported by \code{method="moments"}. Correlation endpoints and zero
variances are allowed; a positive-definite total observation covariance
is still required.

Satterthwaite and asymptotic inference are available. At a whole zero
covariance block, the first-order standard-deviation approximation drops that
block's covariance-uncertainty contribution and its unidentified
correlation; the engine information status reports a boundary
approximation. Kenward--Roger is available for \code{homdiag()} and
homogeneous compound symmetry. It is also available for the equivalent
one- and two-coefficient linear covariance families. Heterogeneous
compound symmetry and AR1 structures with more than two coefficients
have nonlinear covariance restrictions and raise an explicit error if
Kenward--Roger is requested. \code{ddf="Satterthwaite"} and
\code{ddf="asymptotic"} provide inference for those structures. Nonconverged or
unidentified fits do not receive a substituted inference method.

Covariance boundary diagnostics use feasible directions for the chosen
family. A zero block cannot pass merely because derivatives in its
standard-deviation coordinates vanish. Partial zero standard deviations
also require a critical-cone curvature check, including equivalent AR1
correlation signs. If only separated AR1 coefficients have nonzero
variances, a direct check of the first nonzero correlation prevents a
zero derivative of a power of \eqn{\rho} from hiding a decreasing direction.
At zero correlation the check also covers joint changes in correlation and
zero standard deviations, using their leading covariance directions.
These are local necessary optimality checks, not a guarantee of a global
likelihood maximum. If a boundary check cannot certify the requested
structure within its numerical work limits, the fit reports
nonconvergence rather than declaring the boundary optimal. For these
families, \code{object$engine$diagnostics$covariance.boundary.status}
records whether the check was satisfied, a decreasing direction remains,
or a numerical work limit prevented certification.
}
\section{Computation}{Grouping tuples use hash lookups while retaining first-observation
level order. Sparse fixed designs use coefficient-space Woodbury products. Repeated design
rows can also compress the pivoted QR for clearly full-rank designs; ambiguous
rank decisions retain the observation-space QR.
Fitted responses reuse the coefficient BLUPs. Full-rank categorical Type III
tests can use repeated-row compression after checking that the treatment and
sum designs span the same space; other cases retain covariance whitening.}
\section{Covariance optimization}{By default, \code{optimizer="auto"} uses
safeguarded average-information REML (AI-REML) for models with at most 200 free
covariance parameters and only diagonal free covariance blocks. Other models,
including ML fits, use L-BFGS-B. \code{optimizer="AI"} also permits unstructured and the structured random or residual
covariance families described above in REML fits; models above the 200-parameter work limit use L-BFGS-B.
\code{optimizer="L-BFGS-B"} selects L-BFGS-B directly.

Scalar sparse AI-REML models without supplied starts use a feasible residual-only
initial covariance when common-scale profiling is enabled. Supplied parameters,
explicit starts, fixed covariance blocks and zero iteration budgets retain
their requested initialization.

Models with a structured residual family and no supplied covariance starts
or fixed nonzero covariance normalize the common scale of the initial
covariance shape before optimization. The multiplier is \eqn{y'Py/n} for
ML and \eqn{y'Py/(n-r)} for REML. This is a starting-value calculation;
\code{profile.scale} separately controls iterative scale profiling by AI.
Explicit parameter evaluations and \code{maxit=0} retain their requested
coordinates. These residual models allow up to twelve Newton refinement
steps after the main optimizer, using the same objective-decrease and
projected-score acceptance criteria as other models. Refinement does not
relax the convergence tolerance.

AI uses analytic covariance-derivative actions and the fitted projection
operator, without constructing an observation covariance derivative matrix
for every parameter. Each update respects variance bounds and reduces its step
until the likelihood improves. Singular information, failed step reduction or
the \code{ai.maxit} limit transfers the accepted fit to L-BFGS-B with the
remaining main iteration budget. Existing projected-score and positive
semidefinite covariance boundary checks remain in force; a covariance rank
restart uses L-BFGS-B. If the AI continuation fails or fails the score checks,
one additional L-BFGS-B recovery starts from the original feasible parameters
with the same \code{maxit} budget as a fresh L-BFGS-B fit. The better fit is
retained. Reaching the main iteration limit alone does not trigger recovery.
Finite-difference Newton polishing is retained only when
the projected score still exceeds its tolerance. Polishing and
boundary repairs are additional safeguards beyond the main iteration budget.

With \code{profile.scale=TRUE}, AI profiles a common covariance scale whenever
all covariance blocks can be rescaled together and a free residual scale is
available. A nonzero fixed covariance block disables profiling. The conditional
REML scale optimum is \eqn{y'Py/(n-r)} for the current covariance shape, where
\eqn{r} is the estimable fixed-effect rank. Cached factors and operators are
rescaled exactly, scalar sparse scores reuse their trace and quadratic parts,
and the common-scale direction is removed from the AI update.
\code{profile.scale=FALSE} disables this step. This does not assume centered
responses, balanced groups or particular fitted means. Profiling does not
remove an estimated parameter from AIC/BIC or change covariance-parameter
uncertainty calculations. \code{maxit=0} evaluates the supplied starting
coordinates without optimization or profiling.

Engine diagnostics add \code{optimizer}, \code{ai.iterations},
\code{ai.evaluations}, \code{ai.step.reductions}, \code{ai.fallback},
\code{ai.message}, \code{scale.profiled}, \code{scale.profile.updates} and
\code{polish.evaluations}, plus \code{ai.restarts} for the number of recovery
attempts (zero or one). Evaluation counts record covariance evaluations;
scale updates reuse factors and refresh scores separately. These fields
identify fallback and polishing costs. The AI matrix guides optimization only;
Satterthwaite and Kenward--Roger retain their existing information definitions.}
\section{Covariance computation}{The sparse backend supports independent
scalar and correlated random-coefficient terms, diagonal (possibly
heterogeneous) residual covariance, and scalar random terms with relationship
matrices. It forms penalized mixed-model equations in standardized random-effect
space. Minimum-degree ordering and symbolic Cholesky analysis are reused across
optimizer evaluations. Sparse numerical Cholesky and a selected inverse compute
exact likelihoods and scores; no full observation inverse is stored.

For a coefficient covariance \eqn{G=LL'}, the transformed design is \eqn{ZL}.
The symbolic pattern includes all possible nonzero coefficients over the fitted
covariance family. Exactly zero or singular random covariance matrices therefore
retain the same symbolic structure. Covariance derivatives are evaluated with
grouped coefficient-space solves, without inverting \eqn{G} or dividing by its
variances. Temporary derivative work is restricted to batches of grouping
levels. The relationship and REML projection corrections are retained exactly.
This route supports the documented unstructured, diagonal, homogeneous diagonal,
compound-symmetric and first-order autoregressive random-coefficient families,
and fixed positive-semidefinite coefficient covariance matrices.

The sparse route requires the residual covariance to be structurally diagonal
for every permitted covariance parameter value. A correlated coefficient block
with its own supplied relationship matrix, or an observation covariance kernel,
uses an exact covariance backend instead. For designs containing coefficient
blocks, \code{engine$diagnostics$backend.reason} reports the sparse choice or
fallback reason. Forced \code{backend="sparse"} rejects unsupported designs.
Scalar relationship terms may coexist with independent correlated coefficient
terms through the exact relationship correction.

With relationship-only random terms and one positive residual shape, the
whitened relationship Gram matrix is cached for a single kernel as well as
multiple kernels. A complete one-to-one relationship design whose absolute
random-to-residual loading ratio is constant uses its existing relationship
eigenbasis for diagonal variance updates. Every positive eigenvalue is retained;
no ridge, rank approximation or numerical off-diagonal threshold is introduced.
Other designs retain the cached full Gram calculation or the general solver.
Diagnostics \code{relationship.gram.cached} and
\code{relationship.spectral} identify these routes. For the diagonal route,
\code{covariance.factorizations} is zero because variance evaluations use
diagonal weights rather than a new covariance Cholesky factorization.

The block backend detects independent observation components from the
covariance design. It factors those blocks separately and reuses factors only
when the assembled matrices agree exactly after ordering. Ordinary correlated
random-regression and multivariate bases remain supported by this route.

For ordinary multivariate and random-regression terms, the block route also
searches design-based partitions with a small connecting random-coefficient
system. It factors independent observation blocks and applies the exact
Woodbury correction for those connecting effects. Every fitted covariance
direction is retained, including exactly zero or singular random covariance
matrices; no inverse of a random covariance matrix is required. This route
is identified by \code{covariance.backend="block+random"}, with the connecting
coefficient dimension in \code{structured.random.rank}. It currently excludes
supplied relationship and observation kernels. The dense reference remains
available through \code{backend="dense"}.
Independent blocks that agree exactly reuse their factors, and unchanged base
blocks and their design products are retained across variance evaluations.
Covariance-score corrections are contracted in grouped coefficient space;
corrected observation-block inverse matrices are not formed for those scores.

Trait-indicator bases with at most one nonzero per row use that exact support
in covariance and score assembly. General bases retain the full calculation.

The sparse route combines multiple relationship terms in one exact Woodbury
correction, with separate variance scales. For relationship-only models with
one residual scale, it caches the unscaled Gram matrix and uses coefficient-space
solves. Automatic selection keeps the dense reference when multiple relationship
terms have at least as many levels as observations. The block route supports
one scalar relationship term. Singular relationship matrices are supported without a ridge; every
positive eigenmode is retained. Neither balance nor centered responses are
required. Arbitrary kernels and multivariate
relationship terms retain the dense reference route. Automatic selection uses
matrix size and symbolic fill; it is a heuristic, not a timing guarantee.

Verified relationship-kernel product metadata are collected automatically in
\code{model$kernel.products} and \code{engine$diagnostics$kernel.products}.
The detector checks mapped products of other supplied relationship matrices,
then complete ordered Kronecker products. \code{detect.kernel.products=FALSE}
disables detection; \code{kernel.product.tol} sets the relative verification
tolerance (default \code{1e-10}). It preserves the supplied covariance and
selects the computational routes described below when their design requirements
are met. Further details appear in \code{\link{st.kernel.structure}}.

The default \code{backend="auto"} prefers the spectral product solver
for eligible scalar relationship models with one positive diagonal residual
shape. Its default spectral eigen driver prefers DSYEVD. Other model
structures retain the automatic sparse, block or dense fallback.
For a complete crossing with one observation per pair and constant absolute
random-to-residual loading ratio, the solver uses the two factor eigensystems
and tensor transformations. It does not form the full tensor eigenvector matrix.
The tensor shortcut requires reconstruction agreement at floating-point roundoff;
the detection tolerance does not relax this requirement.

For selected crosses, missing responses, replicated observations or nonconstant
loadings, the solver diagonalizes the actual supplied residual-whitened SCA
covariance on observed cells once. It then uses diagonal SCA/residual updates
and an exact lower-dimensional correction for the remaining scalar random terms,
including parental GCA. It never expands a subset to the full crossing and does
not replace a near-product kernel by its reconstruction. Every positive
relationship eigenmode is retained. Zero random variance components and singular
parental covariances are supported without a ridge.

An exactly reconstructed product can also certify PSD through its two parental
factors, avoiding a redundant large PSD decomposition. The factor spectra and
every supplied covariance entry are checked before using this certificate.
Near products retain the full PSD check.

These routes require fewer non-base random-effect levels than observed cells.
The selected spectral route defaults to at most 32766 observed cells, the
largest full DSYEVD order supported by the current 32-bit LAPACK workspace
interface. \code{control$kernel.max.obs} can change this observation limit.
Other designs retain the existing backends. \code{kernel.solver=FALSE} disables computational use while
retaining detection metadata; \code{backend="dense"} selects the reference
solver. Diagnostics \code{kernel.product.solver} report \code{"complete"},
\code{"selected"} or \code{"none"}; \code{kernel.product.correction.rank}
reports the retained correction dimension. The selected route still has a dense
eigendecomposition at setup and is not a matrix-free solver.

\code{control$eigen.solver="dsyevd"} requires DSYEVD explicitly,
\code{"dsyevr"} selects DSYEVR, and \code{"auto"} specifies the default spectral
policy: DSYEVD when its workspace fits the 32-bit interface, and DSYEVR
otherwise. The workspace limit is checked before querying LAPACK. This
fallback can be reached when an explicitly raised observation limit exceeds
32766; an explicit DSYEVD request is rejected at an unsupported order.
Complete crossings retain the tensor shortcut using the factor eigensystems.
The default mapped-product preference uses the user's configured thread count
without changing it. It applies to eligible univariate scalar GCA/SCA models
with three known random-effect covariance matrices, constant coefficient
loadings, one independent residual,
6000--12000 observed responses and at least 16 configured threads. Complete
products retain their spectral shortcut. The rule is based on measured solver
comparisons; larger observation counts retain the prior eligibility rules.
An explicit \code{control$auto.dense.min.threads}, or the option
\code{SelectionTools.auto.dense.min.threads}, overrides this restricted rule
for incomplete products. A finite threshold selects dense when the configured
thread count reaches that threshold; \code{Inf} disables the preference.
Explicit backend choices take precedence.

No additional package or external solver is required. Dense products use R's
BLAS and LAPACK, including the reference implementation. Fixed-effect estimation
uses general GLS. Fixed-effect tests, predictions and their error covariances
use the same fitted covariance operator. Sparse Type III tests use an exact
rectangular whitening. Explicit joint prediction covariance, Kenward--Roger
calculations and parametric bootstrap can still allocate dense intermediate
matrices; sparse fitting does not remove those inference costs.

Scalar covariance blocks in the dense backend use local score accumulators
by default. This reduces result-buffer writes while retaining the same
likelihood and covariance model. \code{control$dense.score.fast=FALSE}
selects the preceding generic score loop for comparison. Blocks with
multiple coefficient dimensions retain that generic loop. The diagnostic
\code{dense.score.fast} records whether the option is enabled.

Diagnostics report \code{covariance.backend}, \code{covariance.blocks},
\code{covariance.max.block}, \code{covariance.factorizations} and
\code{genomic.rank}. For sparse fits, the first two size fields describe the
independent residual rows; the factor count includes the sparse coefficient
factor and any genomic correction. For block fits, it counts distinct
observation-block factors plus any genomic correction; the structured route
reports factors computed in the latest evaluation, excluding reused factors.
Both exclude the
fixed-effect factor. The genomic rank counts retained positive eigenmodes.
\code{sparse.equations} and \code{sparse.factor.nonzeros} report the sparse
coefficient dimension and Cholesky storage, including the diagonal.
\code{sparse.symbolic.analyses} is one per sparse engine context, not per
optimizer evaluation; separate inference reconstructions create new contexts.
These diagnostics do not alter the statistical model.}
\section{Shared printing options}{
\code{st.options()} visibly returns a named list of the current printing
options. Individual values are available through list extraction, for example
\code{st.options()$print.digits}. A saved list can be supplied to
\code{st.options(saved)} to restore the settings. Updates such as
\code{st.options(list(print.digits = 3L))} change only the named entries
and return the complete updated list invisibly. Queries with no argument or
\code{NULL} return the list visibly. The returned list is a snapshot;
changes to that list do not change the settings.

The defaults are \code{print.digits = 4L},
\code{print.char.perline = 75L}, \code{print.top = 10L},
\code{print.bottom = 5L}, \code{print.cov.blocks.top = 3L},
\code{print.cov.blocks.bottom = 1L}, \code{print.ems = FALSE},
\code{print.test.construction = FALSE},
\code{print.tests.for.bare.lsmeans = FALSE}, and
\code{print.intercept.test = FALSE}, \code{print.LSD = TRUE}, and
\code{print.HSD = TRUE}.
Digits range from 1 to 16; the line width ranges from 40 to 10000.
Leading and trailing row and covariance-block counts are nonnegative
integers up to \code{.Machine$integer.max}, or \code{Inf} to display all
rows or blocks in the corresponding section.
The logical switches accept a single nonmissing logical value. Unknown,
duplicate, missing and unnamed option names are errors. Invalid updates
leave every setting unchanged.

SelectionTools printers use these session settings when their corresponding
arguments are omitted or \code{NULL}. Explicit printing arguments take
precedence: \code{digits}, \code{char.per.line}, \code{top} and
\code{bottom} override the four formatting options. For model output,
\code{ems} in \code{...} overrides expected mean square printing for
method-of-moments fits. \code{test.construction} in \code{...} overrides
test construction printing for both likelihood and method-of-moments fits.
The ANOVA type appears where applicable when test construction is displayed.
For least-squares means output, \code{tests}
overrides \code{print.tests.for.bare.lsmeans}; contrasts retain their tests
by default. The switch \code{print.intercept.test} controls stored intercept
tests in ANOVA and fixed-effect test tables; \code{intercept.test} in
\code{...} overrides it for an individual print. Regression coefficient
reports retain the intercept estimate, standard error and test independently.
A fixed-effect test section with no displayed rows is omitted. For an
intercept-only fixed model, it appears when the intercept test is displayed.
The switch does not add an intercept test to moment ANOVA source tables
that contain no such test. Settings affect display only and do not change stored values,
fitting, inference or R's global options. Each newly loaded package DLL
starts with the defaults. The switches \code{print.LSD} and \code{print.HSD}
independently control pairwise difference thresholds beneath least-squares means without
displayed tests. The thresholds remain stored when either switch is false.
}
\section{Method-of-moments output}{Expected mean squares retain full factor
names. Long expressions fill the available line width and continue between
variance components.
Across designs, the ANOVA table contains the numerical results;
optional test-construction blocks identify the numerator and denominator mean
squares. Expected mean squares and test-construction blocks are hidden by
default and are available through the corresponding shared options or
explicit \code{print} arguments. ANOVA and extracted F-test tables omit the
Kind and error-source columns. The optional EMS and test-construction blocks
follow variance or covariance component estimates when these are present.
Trait labels distinguish results for multiple responses, and separate residual
blocks distinguish multiple error strata. Test-construction blocks identify
exact or approximate tests and sources without a requested or available test.
These display choices leave stored ANOVA tables, tests and estimates unchanged.}
\section{Likelihood test construction}{The optional test-construction block
describes the displayed tests after the variance or covariance component
estimates. Exact F tests identify the denominator mean square from the
fitted design when that source was retained. Wald tests identify the
denominator DF method, and likelihood ratio tests identify the comparison
with a reduced model and the reference distribution. The block uses stored
results and requires no additional model fitting.}
\section{Printing extracted results}{The ordinary \code{print} generic uses
the same compiled \code{digits}, \code{top}, \code{bottom} and
\code{char.per.line} controls for models and extracted results. For example,
\code{print(vcov(m), digits=3, top=5, bottom=2, char.per.line=60)} formats the
coefficient covariance matrix. There is no separate printing function.

Numeric results from \code{coef}, \code{vcov}, \code{fitted},
\code{residuals}, \code{model.matrix}, \code{predict}, \code{confint},
\code{logLik}, \code{AIC}, \code{BIC}, \code{nobs} and
\code{df.residual} carry the additional S3 class \code{st.lmm.result}.
\code{AIC} and \code{BIC} return data frames with one row per model, a
\code{DF} column for the parameter count and a column for the criterion.
Formula results retain their formula class and environment. Vectors remain
numeric or integer vectors; matrices retain their dimensions, labels and
matrix inheritance. Numbers retain full precision. The same display class is
used for coefficient and covariance accessors of means and prediction objects.
Existing specialized printers handle their summaries and test tables.

\code{nobs} counts the observed response cells, while \code{df.residual}
subtracts the rank of the full fixed design from this count.
\code{st.lmm.components(object, "observations")} also returns the retained
row count, observed count by trait, original row positions and observed
cell indices in column-major response order. For moments fits,
\code{DF.error.stratum} records the error-stratum DF by trait. These DF
remain the basis for the applicable inference. The observation record
and the model statistics report \code{fitted.type}, which identifies the
target used for the default response residual and its RMSE.
RMSE uses residuals on the original response scale. A weighted moments
fit also reports \code{RMSE.weighted}, based on its precision-weighted
residuals. The \code{fit} matrix in the prediction list has the same layout
as fitted values.

The covariance axes of an intercept-only multivariate block are named
by trait in model order. Covariance blocks with multiple coefficients
within a trait retain their coefficient names. The ordinary residual
block is named \code{units}. Named covariance inputs use the fitted
axis names or the corresponding coefficient names.

New fits also provide this print method for \code{m$engine$sigma},
\code{m$engine$u} and \code{m$engine$diagnostics}, including their extracted
components. It also formats the table returned by \code{anova(m, effects="both")}.
For these lists, \code{top} and \code{bottom} limit displayed components and
the rows within each component; either limit may be \code{Inf}. The model
summary retains its separate three-covariance-block limit.

Printing returns the original object invisibly. It changes no values, labels,
random-number state or global print options. Arithmetic, mathematical functions,
summary operations and ordinary conversions use the original R classes.
Subsetting and transposition retain the result printer. Standard conversions
such as \code{as.numeric}, \code{as.vector}, \code{as.matrix} and
\code{as.data.frame} provide ordinary R values when needed by other code.}
\section{Printed summary}{The ordinary full-model summary begins with the
fixed, random and residual model formulas, followed by convergence and tests.
Fitting controls remain in \code{object$call}. When subset comparison was
requested, the default \code{summary} returns that comparison instead;
\code{summary(object, selection=FALSE)} gives this full-model summary.

When the fixed formula contains numeric predictors, a fixed-effect regression
estimate table follows the random-effect tests and precedes the variance
components. It includes the intercept, if fitted, and terms involving numeric
predictors, with estimates, standard errors, denominator DF, t ratios and
unadjusted two-sided p-values. Purely categorical coefficients are omitted;
\code{st.lsm} provides adjusted factor means. Numeric-by-factor terms are retained:
\code{Diet:Time} can give level-specific slopes, whereas \code{Diet*Time}
contains a reference slope and slope differences with treatment contrasts.
The labels follow the fitted design coding. Reference levels are stated when
the fitted contrast matrix uses a single reference level. Transformations and
matrix-valued numeric bases are included. Multivariate tables are separated by
response. Likelihood fits use the fitted Satterthwaite, Kenward--Roger or
asymptotic method. Ordinary fixed-only moment fits use exact t tests.
Unavailable and non-estimable results are labelled without changing the
inference method. New fits retain this numeric table in \code{object$regression};
it is NULL for models without numeric fixed terms. Older saved fits construct
the table privately when printed, without changing the saved object.

For method-of-moments analyses with a single error term, the ANOVA table
includes a final \code{Residual} row with its degrees of freedom, sum of
squares and mean square. The residual row has no F-test or p-value.
This layout also applies to models with random effects, including nested
villages within regions. Analyses with multiple residual error terms retain
their separate residual table.
The stored tests and fitted estimates are unchanged by this display.

A variance-component table is always printed, including
\code{Residual} (the residual variance) when no random term is specified.
Heterogeneous residual blocks are shown separately. For multivariate responses
and random regressions, row and column labels distinguish variances on the
diagonal from covariances off the diagonal. The session settings
\code{print.cov.blocks.top} and \code{print.cov.blocks.bottom} select the
leading and trailing covariance blocks. An omission count separates the
displayed groups when applicable. A blank line follows each displayed block.
A full labelled covariance matrix is printed when its rows fit
\code{char.per.line} at the requested precision. Wider matrices use a
lower-triangle table. Scalar variances use a separate Component/Variance
table with its own row limits. All estimates are retained in
\code{object$engine$sigma}.
Random-regression coefficients are retained but not printed. The named list
\code{object$engine$u} holds matrices with grouping levels as rows and random
basis coefficients as columns. For example, \code{object$engine$u$Chick}
contains the chick-specific random intercepts and slopes for \code{~Time | Chick}.
These are predicted random deviations (BLUPs). Their sum with the
corresponding fixed regression contribution gives an individual regression. A random-slope-only
model has no random-intercept column.

Fits criteria and model statistics are combined at the end in two side-by-side
columns of labelled values. Multivariate summaries show the overall fit criteria
once and identify each trait's statistics separately. Numbers use fixed decimal
notation; small p-values and large F-statistics are displayed as bounds at the
stated thresholds. Random-test sections are omitted when tests were not requested
or there are no random terms; variance estimates are always shown. Printing does
not modify the underlying estimates or tests. Printed headings use \code{DF1} and \code{DF2} for numerator and denominator degrees of freedom, \code{DF} for a single degrees-of-freedom column, and \code{F}, \code{t}, and \code{p} for test statistics and probabilities. Wrapped model formulas use hanging indentation; formulas in a random-term list are shown on separate lines. The data frames returned by \code{anova} use the same compiled printer while retaining their stored column names.}
\section{Designed experiments and denominator DF}{For balanced split plots,
Satterthwaite and Kenward--Roger reproduce the exact error-stratum DF when
the corresponding variance components are estimable. In the balanced split-plot
datasets \code{\link{v-lmm-alfalfa}} and \code{\link{v-lmm-oats}},
the variety effect has 10 denominator DF, and subplot effects and their
interactions have 45. An explicitly requested \code{ddf="asymptotic"} instead
uses an infinite denominator DF.

SAS containment uses rank contributions of syntactically containing random terms,
or \eqn{n-\mathrm{rank}[X\ Z]} if none contains the fixed term. It is not
currently an \code{st.lmm} option. In incomplete-block and lattice mixed models,
Satterthwaite and Kenward--Roger can have fractional DF that differ from containment
or a least-squares analysis treating blocks as fixed. Comparable results require
the same model and DF method.}
\section{External inference interfaces}{
When users install \pkg{emmeans} or \pkg{multcomp}, their inference functions
can consume supported \code{st.lmm} objects directly. SelectionTools does not
load these packages and does not require them for installation or fitting.
Methods become available when the corresponding external namespace is loaded.
All SelectionTools-side calculations use compiled C routines with thin R entry
points. No interface refits the model or calls an optimizer.

For univariate likelihood fits, \code{emmeans::emmeans(object, ~ factor)}
uses retained predictor data, the compiled fixed-effect design, covariance and
estimability information. The optional \code{ddf} argument selects
\code{"Satterthwaite"}, \code{"Kenward-Roger"}, or \code{"asymptotic"}; the
fitted object's method is the default. Degrees of freedom are calculated for
each requested linear combination. Kenward--Roger requires REML and uses the
matching adjusted covariance. Rank-deficient designs retain a non-estimability
basis. Unavailable finite-sample inference is not replaced by asymptotic inference.

\code{multcomp::glht(object, linfct = K)} and supported
\code{multcomp::mcp} comparisons use fixed coefficients and their matching
covariance. By default \code{glht} uses its normal approximation. A user may
supply one finite nonnegative \code{df} at the \pkg{multcomp} interface; positive values request its common-DF
t inference, and zero requests the normal approximation. Contrast-specific
Satterthwaite degrees of freedom are not collapsed into a common value.

The first implementation excludes multivariate fits, moments fits, selected-model
collections, emmeans data/covariance overrides, argument-supplied offsets, and
glht coefficient/covariance overrides. Formula offsets are supported by emmeans.
The glht interface rejects rank-deficient fits; the estimability-aware
emmeans interface supports those fits. Interaction comparisons follow the external
package's own contrast and averaging conventions. Automatic trends, transformations,
and advanced external options require validation for the model concerned.
}

\value{\code{st.options} returns a named list of the ten current printing
options: visibly for queries with no argument or \code{NULL}, and invisibly
after an update supplied as a list. The model-fitting function returns
a serializable S3 object of class \code{st.lmm}, with the call/formulas/data, compiled \code{engine} fit, compiled \code{model} metadata/designs, cached \code{tests}, numeric fixed-regression \code{regression} table, per-trait \code{statistics}, and inference controls. The engine includes coefficient covariance, covariance blocks, random-effect predictions, alias/estimability maps and parameter-information derivatives for subsequent linear-function inference.
When subset comparison is requested, the object also contains
\code{selection}, holding source-column groups, retained models, comparison
statistics, search settings and diagnostics. Its default summary has class
\code{st.lmm.selection.summary}. The full-model fit is not replaced by a
selected candidate.}
\seealso{\code{\link{st.lmm.profile}}, \code{\link{st.predict}}, \code{\link{st.lsm}}, \code{\link{st.get.A}}}

\examples{
saved.options <- st.options()
st.options(list(print.digits = 3L, print.ems = TRUE))
st.options(saved.options)
}
