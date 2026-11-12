\name{st.lmm.by}
\alias{st.lmm.by}
\alias{summary.st.lmm.by}
\alias{coef.st.lmm.by}
\alias{confint.st.lmm.by}
\alias{print.st.lmm.by}
\alias{print.summary.st.lmm.by}
\title{Fits Independent Gaussian Models to Groups of Observations}
\description{
Observations are split by a grouping variable, and the same model specification is fitted
separately in each group.  The result retains each fitted model, its original row mapping,
warnings and fitting status.  Fixed-only univariate regressions can optionally
use a common pooled residual variance for coefficient inference.
}
\usage{
st.lmm.by(fixed, group, data, random = NULL, cov = ~units,
          method = "REML", ddf = "Satterthwaite", type = 3L,
          random.test = "none", control = list(), ucov = NULL,
          weights = NULL, contrasts = NULL, subset = NULL,
          na.action = NULL, offset = NULL, pool = FALSE)

\method{summary}{st.lmm.by}(object, pool = NULL, level = .95, ...)
\method{coef}{st.lmm.by}(object, ...)
\method{confint}{st.lmm.by}(object, parm, level = .95, pool = NULL, ...)
\method{print}{st.lmm.by}(x, digits = NULL, top = NULL, bottom = NULL,
                         char.per.line = NULL, ...)
\method{print}{summary.st.lmm.by}(x, digits = NULL, top = NULL,
                                 bottom = NULL, char.per.line = NULL, ...)
}
\arguments{
  \item{fixed}{A two-sided fixed-effect formula, as in \code{\link{st.lmm}}.}
  \item{group}{A single data-column name, a formula of the form
  \code{~column}, or an atomic grouping vector with one value per original
  data row.  A column symbol can be used directly.  Groups are ordered by
  their first appearance after subsetting.  Selected group values must not
  be missing or nonfinite.}
  \item{data}{An R data frame or other data representation accepted by
  \code{\link{st.lmm}}, including a bitpacked SelectionTools marker data
  set or a SelectionTools hybrid data set.}
  \item{random, cov, ddf, type, random.test, control, ucov, contrasts}{Model
  specification and fitting arguments passed to each independent fit.
  Further details appear in \code{\link{st.lmm}}.}
  \item{method}{\code{"ML"} or \code{"REML"}.  Method-of-moments fitting
  is not supported by this grouped interface.}
  \item{weights}{Optional numeric observation precision weights, evaluated
  on the original data and then sliced for each group.  The validation and
  zero-weight rules of \code{\link{st.lmm}} apply within each group.}
  \item{subset}{Optional logical or integer row selection, evaluated on the
  original data before grouping.  Repeated integer indices are retained.}
  \item{na.action}{A standard missing-value action, applied separately within
  each group.  \code{NULL} retains the compatibility behavior of
  \code{\link{st.lmm}}.}
  \item{offset}{An optional known numeric offset, evaluated on the original
  data before grouping.  A scalar is recycled.  It is added to any offsets
  specified in the fixed-effect formula.}
  \item{pool}{Logical for \code{st.lmm.by}.  TRUE requests coefficient inference
  using a common residual variance across groups.  This requires fixed-only
  univariate models with \code{cov = ~units}.  In \code{summary} and
  \code{confint}, \code{NULL} uses the choice made when fitting; \code{TRUE}
  or \code{FALSE} requests pooled or separate-model inference.}




  \item{object, x}{An \code{st.lmm.by} object, or its summary for the
  corresponding print method.}
  \item{parm}{Optional coefficient names or indices into the columns of
  \code{coef(object)}.  Omission selects all coefficients.}
  \item{level}{The confidence level, strictly between zero and one.}
  \item{digits, top, bottom, char.per.line}{Per-call display controls, as in
  \code{\link{st.lmm}}.}
  \item{...}{Additional generic-method arguments; currently unused.}
}
\details{
Repeated fits run in the current R process. Numerical kernels use the
existing OpenMP and BLAS thread controls; \code{\link{st.set.num.threads}} describes these controls.
R callbacks execute on the main R thread. Intermediate state stays in memory.

The models are independent fits.  In particular, a random-effect level with
the same label in two groups does not share an estimated random coefficient
between those groups.  Fixed-effect estimates and covariance parameters are
estimated separately.  Shared effects across groups are represented by the
grouping structure in a single \code{\link{st.lmm}} model.

Grouping, weights, the argument offset and the subset are evaluated on the
original data.  The subset is applied before splitting.  Each subgroup is
then treated as a fresh data set: transformations such as \code{scale(x)}
or \code{poly(x, 2)} are learned from that subgroup, before its missing-value
action.  Declared factor levels and contrast attributes are retained as in
ordinary data-frame subsetting.  Aliased coefficients remain explicitly
marked as not estimable.  The prepared subgroup designs, evaluated inputs,
contrasts and referenced formula bindings are retained, so saved tasks do
not depend on later changes to the original data or caller bindings.

Marker or hybrid data and any requested relationship matrices are resolved
before fitting the groups.  The retained subgroup requests can therefore
be evaluated after the corresponding in-memory data set has been removed.
The original marker or hybrid representation remains unchanged.

An error in an individual group's preparation or fit is retained in that
group's result.  A fit that does not converge is retained with a separate
status.  Groups are never silently dropped, and coefficient inference is
unavailable for unsuccessful fits.  The result's \code{complete} field is
true only when every group fit converged.  Each ordinary fitted model uses
its local data-row mapping; the grouped result additionally records original
data-row indices and positions in the selected layout, including duplicates.

Pooling assumes a common residual variance after accounting for any known
precision weights.  It is available only for fixed-only univariate models
with one independent residual-scale component. Here \eqn{e_g} is the
residual vector in group \eqn{g}, \eqn{W_g} is its diagonal precision matrix,
\eqn{n_g} is its number of observations and \eqn{r_g} is its fixed-design rank.
The pooled variance is
\deqn{\widehat{\sigma}^2 =
  \frac{\sum_g e_g^\mathsf{T} W_g e_g}{\sum_g(n_g-r_g)}.}
The coefficient covariance in each group is this variance times the inverse
weighted fixed-design information on its estimable columns.  Confidence
intervals and tests use the Student distribution with
\eqn{\sum_g(n_g-r_g)} degrees of freedom.  No pooled inference is produced
unless all required fits succeed and the pooled residual degrees of freedom
are positive.

\code{summary} adds a coefficient table with estimates, standard errors,
degrees of freedom, test statistics, probabilities, confidence limits and
inference statuses.  It retains the group inventory and reports whether
pooled inference was requested.  \code{coef} returns a matrix with one row
per group and the union of coefficient names as columns; unavailable or
aliased coefficients are \code{NA}.  \code{confint} returns the coefficient
table restricted to the selected coefficients, retaining explicit rows for
failed groups.  The interval level and inference interpretation are stored
as attributes.  These methods do not modify the fitted object.
}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lmm.refit}},
\code{\link{st.lmm.components}}}
