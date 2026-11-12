\name{st.lmm.components}
\alias{st.lmm.components}
\alias{st.lmm.singular}
\alias{st.lmm.pca}
\title{Components and Random-Coefficient Covariance Diagnostics of a Fitted Model}
\description{
Extracts coefficients, covariance components, model construction and
estimability information from a Gaussian linear mixed model. The related
functions diagnose the numerical rank of each random-coefficient covariance
block and return its principal components. Construction and calculations are
performed in C.
}
\usage{
st.lmm.components(object, name = NULL, traits = NULL, block = NULL,
  term = NULL)
st.lmm.singular(object, tol = 1e-4)
st.lmm.pca(object, tol = 1e-7)
}
\arguments{
\item{object}{A univariate or multivariate \code{st.lmm} fit.
\code{st.lmm.components} accepts ML, REML and moments fits.
\code{st.lmm.singular} and \code{st.lmm.pca} accept ML and REML fits.}
\item{name}{NULL for the standard component list, or a character vector of
component names listed below. A single name returns the component itself.
Several names return a named list in the requested order. Names must be
nonempty, distinct and nonmissing. Matching is exact.}
\item{traits}{NULL or a character vector of distinct fitted response trait
names. The variance, covariance and correlation components support this
selection. An explicit selection uses the supplied trait order; NULL
uses fitted order.}
\item{block}{NULL for all covariance blocks, or one exact fitted covariance
block name. The variance, covariance and correlation components support
this selection. The result for one selected covariance or correlation
block is its matrix. The variance matrix retains one row for the block.}
\item{term}{One exact, nonmissing and nonempty hypothesis term name for
\code{sscp.H}. This argument is required for \code{sscp.H}, including
requests with several component names. It is NULL for other requests.}
\item{tol}{A finite number strictly between zero and one. A principal
component is counted in the numerical rank when its standard deviation,
divided by the square root of the reported reference variance, is strictly
greater than \code{tol}. Further details appear in Details for the reference scale.}
}
\details{
\subsection{Components}{
The standard result contains the following entries.
\describe{
\item{\code{fixed}}{Named fixed-coefficient estimates in the full original
model-matrix order. Aliased columns retain their \code{NA} placeholders.
This is also the default coefficient content of \code{coef(object)}.}
\item{\code{random}}{A named list of random-coefficient predictions. An
ordinary random block is a matrix with group levels in rows and coefficient
names in columns. A direct covariance-kernel block has a single column with
one prediction per response cell. Moments fits and fits without random terms
return an empty list.}
\item{\code{vcov.fixed}}{The model-based fixed-coefficient covariance matrix
at the fitted covariance parameters. It is unadjusted even when the fit
requested Kenward--Roger inference. This extractor does not change the
adjustment policy used by \code{vcov(object)}.}
\item{\code{covariance}}{A list with \code{random} and \code{residual}
entries, each containing the fitted coefficient covariance matrices in
their original coefficient units. The axes of an intercept-only
multivariate block identify the traits. In a likelihood fit, known covariance
between grouping levels contributes to the full \code{G} matrix. A moments
fit retains the unconstrained covariance component estimates. With one
\code{block}, the component contains the selected matrix directly.
Trait selection retains every coefficient of
each selected trait, in coefficient order within the trait. The grouped
result omits blocks without a selected coefficient. A single selected block
requires at least one matching coefficient.}
\item{\code{residual.scale}}{The square root of the residual variance
multiplier when there is one univariate, one-dimensional residual covariance
block. Otherwise \code{NA_real_}. With precision weights or a residual basis,
the multiplier need not equal every observation's residual standard
deviation; \code{R} gives the observation-level residual covariance.}
\item{\code{ngroups}}{Named counts of the levels represented in each random
covariance block, including extra levels supplied through a known covariance
matrix. The \code{unit} attribute distinguishes group levels from response
cells in direct covariance-kernel blocks.}
\item{\code{joint}}{The stored joint fixed/random prediction-error covariance
matrix, or a zero-by-zero matrix when it was not stored. Explicitly requesting
\code{name="joint"} reconstructs this matrix at the fitted covariance
parameters when necessary for ML and REML fits. Reconstruction holds
the fitted parameters fixed and leaves \code{object} unchanged. In a moments
fit this component is a zero-by-zero matrix.}
\item{\code{joint.order}}{A data frame declaring every position in
\code{joint}: \code{position}, a unique \code{label}, \code{kind}
(fixed or random), \code{term}, \code{level}, \code{coefficient},
\code{fixed.column}, \code{block} and \code{group.index}. Inapplicable
entries are \code{NA}. The fixed positions are independent columns in the
engine's pivot order. In a moments fit this data frame has zero rows.}
\item{\code{estimable}}{A list with logical \code{fixed} flags for individual
coefficients in the original parameterization, the independent
\code{fixed.columns}, the \code{fixed.map} relating original columns to
the independent design, and logical \code{fitted} flags for response-cell
means. A coefficient can have a finite representative estimate but fail to
be individually estimable in an aliased parameterization. Fitted flags are
\code{NA} for rows omitted because model predictors are unavailable. The
\code{scope} entry explicitly separates these estimability flags from fit
convergence and covariance-parameter identification.}
\item{\code{metadata}}{The method, formulas, response traits, observed
response-cell indices, omission mask, offsets, weights, fixed rank and
aliased column names, convergence state and engine diagnostics, available
covariance-parameter information rank and status, joint-covariance storage
status and interpretation, matrix ordering, and fitted-object version.
The \code{availability} entry describes the status of random predictions,
joint covariance, its ordering and \code{G}.
Unavailable covariance-parameter information is reported as such. Information
rank deficiency is distinct from fixed-column aliasing and from convergence
failure. It can also reflect a boundary parameterization, so
\code{covariance.parameters$confounded} is \code{NA}: these extractors do
not separately diagnose covariance-parameter confounding.}
}
Additional named components are available:
\describe{
\item{\code{series}}{The second-stage inputs and target declaration for a
series fit. The \code{data} data frame and \code{rcov} matrix are the
retained second-stage inputs. Their row names and any covariance-axis labels
are preserved. The matrix is the inherited or supplied target covariance on
the response scale. The \code{cells} data frame has the complete target
labels, trial, genotype, trait and data-row mapping before second-stage
row selection. The \code{reference} entry contains the declared grids and
weights, or NULL for automatic targets. The \code{genotype} entry identifies
the resolved genotype column. The \code{covariance.provenance} entry records
an explicitly supplied joint target covariance and is NULL under the
default independence assumption. The \code{independence} entry records the
earlier fits, observation counts and disjointness check. A fit without a
second-stage series declaration gives an error.}
\item{\code{variance}}{A matrix with one row per covariance block and one
column per trait, including when a single block or trait is selected.
The entries are diagonal elements of the fitted coefficient covariance
matrices. A trait with no coefficient in a block has an \code{NA}
entry. A selected trait with several coefficients in a block requires
the \code{covariance} component to obtain their separate variances.
For a \code{homdiag} block with a common coefficient variance, the
matrix gives that shared variance.}
\item{\code{correlation}}{The correlation matrices of the fitted
coefficient covariance blocks, with the same structure and selection
rules as \code{covariance}. A coefficient with a nonpositive or nonfinite
variance has \code{NA} correlations. Positive variances are standardized
using the fitted covariances, including unconstrained moments estimates.}
\item{\code{observations}}{Counts of retained data rows and observed
response cells, overall and by trait, the fixed-design rank and residual DF,
and the row and cell indices. Moments fits also provide the residual error
stratum DF by trait.}
\item{\code{tests}}{A table of the stored fixed, random and joint tests in
a common format, including their methods, statistics, DF and p values.
The joint rows of \code{anova(object, effects = "all")} include the
MANOVA statistics and their F reference calculations.}
\item{\code{sscp}}{A list of sums-of-squares-and-products results for
method-of-moments fits. \code{H} is the list of MANOVA hypothesis matrices
named by model term, and \code{E} is the error matrix for the original
fit. \code{strata} and \code{pooled} contain the original and pooled
orthogonal source matrices, respectively, named by source. Matrices have
traits on both axes, including one-by-one matrices for a single trait.
An unavailable matrix has zero rows and columns, and an unavailable list
has zero elements.}
\item{\code{sscp.H}}{The stored hypothesis sums-of-squares-and-products
matrix for \code{term}, with the trait names on both axes. A term without
a stored hypothesis matrix gives an error.}
\item{\code{sscp.E}}{The error sums-of-squares-and-products matrix for the
original fit, with the trait names on both axes. A univariate error
matrix is one by one. A fit without an error matrix gives a
zero-by-zero matrix.}
\item{\code{predictions}}{The stored hybrid prediction data frame with two
parental identifier columns named for the fitted parental pools and an
integer \code{source.row} column. For \code{prediction = "NAonly"},
\code{source.row} identifies the row in the retained performance data.
For \code{prediction = "Full"}, the parental identifiers specify the
cross and \code{source.row} is \code{NA}. The table has a numeric
\code{prediction} column for one trait, or a \code{prediction.<trait>}
column for each fitted trait. For \code{"NAonly"}, an observed trait cell
has no cached prediction and its entry is \code{NA}. Extraction reuses the
saved predictions. A non-hybrid fit returns a zero-row data frame with
factor columns \code{parent1} and \code{parent2}, integer \code{source.row}
and numeric \code{prediction}.}
}
The following dense matrices are returned only when explicitly named:
\describe{
\item{\code{X}}{The complete fixed-effects design, including aliased columns.}
\item{\code{Z}}{The combined random-coefficient design.}
\item{\code{G}}{The full covariance of the random coefficients. For an
ordinary block with coefficient covariance \eqn{S} and known level covariance
\eqn{A}, its block is \eqn{S\otimes A}; independent levels use the identity
matrix for \eqn{A}. A direct kernel block uses an identity design in
\code{Z} and its fitted kernel covariance in \code{G}. Moments fits with
random terms return a zero-by-zero matrix for \code{G}. Their component
estimates are available in \code{covariance}.}
\item{\code{R}}{The residual covariance between response cells, including
all residual blocks, their bases, any level covariances, and precision
weights.}
}
These matrices can be large. The standard component list does not construct
\code{X}, \code{Z}, \code{G} or \code{R}; explicitly requesting an uncached
\code{joint} matrix also requires dense storage for that matrix.
}
\subsection{Ordering, units and uncertainty}{
The rows of \code{X}, \code{Z} and \code{R} follow the retained model's
response-cell order, with all rows of the first trait before all rows of the
next trait for a multivariate response. The \code{observed} indices in
\code{metadata} identify the response cells used in fitting. Matrix rows
can therefore include cells without an observed response. Offsets are
reported separately in \code{metadata$offset}; they are not columns of
\code{X}.

Within a random covariance block, coefficient positions vary first over all
group levels of the first coefficient, then all levels of the next
coefficient. Thus they match the column-major vectorization of the
corresponding matrix in \code{random}. Blocks remain separate even if they
refer to the same grouping variable. \code{Z}, \code{G} and the random
positions of \code{joint} share this order. Matrix names and
\code{joint.order} declare it explicitly.

With covariance parameters held at their fitted values, \code{joint}
estimates
\deqn{\operatorname{Cov}\left[\begin{pmatrix}\widehat\beta-\beta\\
\widehat u-u\end{pmatrix}\right].}
The random block is prediction-error covariance with fixed effects estimated.
It includes the contribution of fixed-effect estimation and the fixed/random
cross covariance. It is not the conditional covariance calculated with
fixed coefficients treated as known, nor does it include covariance-parameter
uncertainty or residual variation of a future observation. Only independent
fixed columns occur in \code{joint}; \code{estimable} supports interpretation of
an aliased fixed model. The returned matrix carries
\code{uncertainty.kind} and
\code{includes.covariance.parameter.uncertainty=FALSE} attributes.
}
\subsection{Numerical rank and principal components}{
For each random-coefficient covariance matrix \eqn{S_k}, the principal
component calculation returns its eigenvalues in descending order and an
orthonormal eigenvector matrix. Both covariance and eigenvalues use the
original coefficient units. This calculation concerns \eqn{S_k}, not the
rank of a supplied level covariance matrix \eqn{A_k}. A singular supplied
relationship matrix alone does not make \code{st.lmm.singular} return TRUE.

The quantity \eqn{\bar r} is the mean diagonal of \eqn{R} over observed
response cells, and \eqn{\lambda_{k,1}} is the largest eigenvalue of \eqn{S_k}.
The reported reference variance is
\eqn{s_k=\max(\bar r,\lambda_{k,1})}. The numerical rank is the number of
eigenvalues satisfying \eqn{\sqrt{\lambda_{k,j}/s_k}>\mathrm{tol}}.
If the reference variance is zero, every relative standard deviation is
defined to be zero. Negative eigenvalues within numerical roundoff are
truncated to zero; a materially indefinite covariance matrix is an error.

This definition is unchanged by a common rescaling of the responses and is
specified in the original random-coefficient basis. Changing predictor units
can change a relative-rank assessment. The definition differs from a
boundary check on optimizer or Cholesky coordinates. The eigenvalues,
reference variances, tolerances and ranks are returned so the criterion can
be inspected directly.

An eigenvector's sign is arbitrary. Repeated eigenvalues identify a principal
subspace rather than unique eigenvectors. Comparisons should use eigenvalues
and subspaces instead of requiring identical eigenvector signs or bases.
Fixed-column aliasing, covariance-parameter information deficiency and lack
of convergence are reported separately and do not themselves set the
singularity flag. The convergence attribute should be considered when
interpreting a diagnostic from an unfinished fit.
}
}
\value{
\code{st.lmm.components} returns the value of a single named component.
Multiple names or a NULL request return a named list of components.
Matrices retain their dimensions, including one-row and one-column cases.
Extracted values retain the package's result-display class in addition to
their ordinary R structure.

\code{st.lmm.pca} returns a named list with one entry per random covariance
block. Each entry contains \code{sdev}, \code{rotation},
\code{eigenvalues}, \code{relative.sdev}, \code{rank}, \code{dimension},
\code{reference.variance}, \code{tol} and \code{singular}. The list has
\code{tol}, \code{rank.definition}, \code{mean.residual.variance} and
\code{converged} attributes. A model without random blocks returns an empty
list with these attributes.

\code{st.lmm.singular} returns one logical value: TRUE if any random
coefficient covariance block has numerical rank smaller than its dimension.
Named attributes \code{terms}, \code{rank}, \code{dimension},
\code{eigenvalues}, \code{relative.sdev} and \code{reference.variance}
give the block-level results. It also carries \code{tol},
\code{rank.definition}, \code{mean.residual.variance} and
\code{converged}. A model without random blocks returns FALSE.
}
\seealso{\code{\link{st.lmm}}, \code{\link{st.predict}},
\code{\link{st.lsm}}}
\keyword{models}
