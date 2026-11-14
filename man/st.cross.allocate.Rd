\name{st.cross.allocate}
\alias{st.cross.allocate}
\alias{print.st.cross.allocation}
\title{Integer Progeny Allocation to Fixed Crosses}
\description{
Allocates a fixed produced-progeny budget across a fixed set of unordered
parent pairs, maximizing predicted mean after within-family selection.
Retained counts are fixed; produced counts are integer decision variables.
}
\usage{
st.cross.allocate(plan, scores, total.progeny, n.selected = 10,
                  min.progeny = n.selected, max.progeny = total.progeny,
                  accuracy = 1, control = list(),
                  selection.method = c("burrows", "order.statistics"),
                  constraints = NULL, genetic.model = NULL)

\method{print}{st.cross.allocation}(x, digits = NULL, top = NULL,
                                  bottom = NULL, char.per.line = NULL, ...)
}
\arguments{
  \item{plan}{An \code{st.cross.plan} result, or a data frame containing
  nonempty character columns \code{P1Name} and \code{P2Name}. Self crosses
  are excluded. Pair order and orientation follow their first input
  occurrence. Repeated unordered pairs are pooled as described below.}
  \item{scores}{A data frame with character \code{P1Name}, \code{P2Name}
  and finite numeric \code{mu}, \code{va} columns. Variances must be
  nonnegative. Every pooled planned pair has exactly one score row;
  additional unplanned score rows are permitted. Existing \code{es}
  values are not used because selection intensity depends on the new size.}
  \item{total.progeny}{Positive supported whole-number produced budget.}
  \item{n.selected}{Positive retained counts, supplied as a numeric scalar
  or one value per input plan row.}
  \item{min.progeny, max.progeny}{Positive integer produced-count bounds,
  each a numeric scalar or one value per input row. Every row satisfies
  \code{n.selected <= min.progeny <= max.progeny}.}
  \item{accuracy}{Finite prediction accuracy in \eqn{[0,1]}, a numeric
  scalar or one value per input row. Repeated unordered pairs must have
  identical accuracies.}
  \item{selection.method}{Selection-response calculation. The default
  \code{"burrows"} retains the finite-family approximation used by
  \code{gs.cross.eval.es}. The opt-in \code{"order.statistics"} method
  numerically evaluates the expected mean of the retained top normal
  order statistics for the actual produced and retained counts.}
  \item{constraints}{Optional practical restrictions described in
  \code{\link{st.cross.constraints}}. Explicit restrictions intersect
  inherited ones. Produced-count restrictions can couple families and
  require constrained search; the separable allocation optimum then does
  not certify the constrained problem. Directions, when needed, are
  identified explicitly.}
  \item{control}{Named list with unique supported components:
  \describe{
    \item{\code{max.cells}}{Positive integer limit on the estimated DP
    utility, traceback and numeric-row cells; default \code{10000000}.
    The order-statistic method estimates its quadrature, cache and
    allocator workspaces as numeric-equivalent cells. This is a search
    workspace estimate, not a total byte-memory limit; mandatory output
    evaluation can still require work after a search limit is reached.}
    \item{\code{max.transitions}}{Positive integer limit on a conservative
    DP transition-count estimate; default \code{100000000}. For the
    order-statistic greedy path, the cap counts allocated extra
    individuals.}
    \item{\code{tolerance}}{Nonnegative response-gap tolerance; default
    \code{1e-10}. A completed search is qualified optimal when its
    enclosed response gap is at most this value times
    \code{max(1, abs(response))}.}
    \item{\code{algorithm}}{For \code{selection.method="order.statistics"},
    one of \code{"auto"}, \code{"greedy"} or \code{"dp"}. The automatic
    path uses discrete marginal gains; the dynamic program uses the same
    order-statistic objective. This component is not a Burrows
    control.}
  }}
  \item{genetic.model}{Optional additive DH genetic model, as described in
  \code{\link{st.cross.evaluate}}. The model supplies approximate
  selection-aware founder contributions and optional final actual-genotype
  scenarios. An explicit \code{future} group also requests the final
  next-cycle assessment described in \code{\link{st.cross.future}}, without
  changing this allocator's current-generation objective.
  A source core and cap are inherited with \code{NULL}; scenario
  execution always requires a current request. An active
  \code{selected.founder.theta} requires an original supplied or inherited
  \code{K} and couples the actual produced and retained counts.}
  \item{x}{An \code{st.cross.allocation} result.}
  \item{digits, top, bottom, char.per.line}{Print controls. \code{NULL}
  uses the corresponding \code{st.options()} setting. Top and bottom
  counts may be \code{Inf}.}
  \item{...}{Reserved for S3 printing.}
}
\details{
For pooled family \eqn{e}, \eqn{r_e} is the fixed retained count,
\eqn{g_e} the produced count, \eqn{a_e} the accuracy, and
\eqn{R=\sum_e r_e}. The predicted selected mean is
\deqn{m_e(g_e)=\mu_e+a_e\sqrt{v_e}\,i_B(r_e,g_e),}
and the objective is \eqn{\sum_e(r_e/R)m_e(g_e)} subject to integer
family bounds and a produced budget \eqn{\sum_e g_e=B}.
The default uses the same finite-family Burrows selection-intensity
calculation as \code{gs.cross.eval.es(N=r_e, G=g_e)}. It uses zero intensity
when all produced progeny are retained. Accuracy scales the selection
response rather than the family mean.

For \code{selection.method="order.statistics"}, the intensity is
\deqn{i_O(r,g)=\frac{1}{r}\sum_{j=g-r+1}^{g}E[Z_{(j:g)}],}
where \eqn{Z_{(j:g)}} is a standard-normal order statistic. This replaces
\eqn{i_B} in the selected-mean formula. The finite-normal expectation
follows Okoye, Frisch and Herzog (2026), equations (4)--(8). Retained-count
weights also cover unequal fixed retained counts; they reduce to the
paper's unweighted family average when retained counts are equal.
The expectation is exact under the normal scoring model; its numerical
evaluation and the fit of that model to biological progeny are separate
qualifications. For phenotypic selection, \code{accuracy} is the square
root of heritability under the within-family regression model.

The order-statistic greedy path compares actual one-individual response
increments, rather than derivatives of a continuous family-size extension.
These increments decrease with produced count for a fixed retained count.
Consequently, choosing the largest available weighted increment solves
the separable integer problem in exact arithmetic. Bounds truncate the
increment sequences. This argument does not extend to joint pair choice,
changing retained counts or constraints that couple produced family sizes.
The alternative dynamic program evaluates the same requested score.

Repeated reciprocal or identically oriented pair rows describe a single
biological selection family. Their retained counts and produced lower and
upper bounds are added before evaluating selection intensity. Selection
therefore retains the pooled total from the pooled family. It does not
perform independent selection in the original repeated rows.

For the default Burrows method, a two-row integer budget dynamic program
enumerates all feasible family
sizes; integer traceback recovers the allocation. The constant weighted
family-mean component is excluded from DP comparisons. No concavity or
continuous-size assumption is required. The algorithm has time complexity
\eqn{O(FBW)} and traceback space \eqn{O(FB)}, where \eqn{F} is the number
of pooled families, \eqn{B} the budget above summed minima, and \eqn{W}
the largest feasible size range. These are fixed-pair allocation bounds,
not claims of a joint optimum over parent pairs and progeny sizes.

When \code{plan} is an equal-family \code{st.cross.plan} result without
explicit contribution-basis provenance, nominal retained
parental contributions must satisfy its inherited contribution bounds,
fixed-contribution requirements and relatedness limit. A retained family
contributes half its nominal gene share to each parent. The inherited
\code{K} is used to report retained and produced group relatedness
separately. Produced contributions are descriptive and are not constrained
by retained-share limits. A family-slot diversity constraint inherited
from the original fixed-pair plan is unchanged by altering produced sizes.
A plain data frame has no inherited \code{K} or planning constraints, so
its relatedness values are \code{NA}.

A source plan from the optional \code{progeny} model explicitly records
\code{family.slots} as the basis of its top-level parental and relatedness
restrictions. Those inherited restrictions retain that basis, even when
this allocation supplies unequal retained counts. Reported retained and
produced contributions and relatedness still use the actual corresponding
counts. Explicit practical restrictions retain their own declared bases.
The allocator's \code{total.progeny}, retained counts and produced bounds
define a new fixed-pair domain; they replace the source plan's production
design, including source equal-size bounds. Its scoring method follows this
call's \code{selection.method}. Original parent, allowed-pair, reproductive
role and practical restrictions remain in force.

Neither nominal contribution vector estimates realized genomic ancestry
after within-family selection. Crossing populations with the returned
variable sizes uses the existing low-level crossing commands; the
equal-family execution adapter is unchanged.

A completed DP also propagates an upward addition enclosure over the
computed binary64 response table. The reconstructed response receives a
downward addition enclosure. Their difference is the reported response gap.
This qualification concerns arithmetic on evaluated utilities and does
not bound error of normal-distribution or finite-family approximations.
Status \code{optimal} requires completion and the requested gap tolerance;
otherwise a completed search has status \code{numerically.unresolved}.
If either resource estimate exceeds its cap, a deterministic feasible
allocation fills family capacity in input order and is returned with
status \code{resource.limit}, infinite objective bound and infinite gap.
Such a result is an incumbent without an optimality claim.

The order-statistic method additionally reports estimated numerical
integration uncertainty and identifies its scoring method and allocator.
Status \code{numerically.optimal} means that the estimated numerical
response gap meets the requested tolerance. Quadrature error estimates
are not rigorous interval enclosures. A larger unresolved estimate is
reported as \code{numerically.unresolved}; a resource-limited result is
reported as \code{resource.limit}. These qualifications do not certify
the finite-normal model for arbitrary biological progeny distributions.
A limited greedy search retains the allocated increments and fills the
remaining feasible budget in input order. A workspace-limited search
can return the initial feasible allocation before greedy increments.
The optional genetic model and its inherited selected-founder cap use actual
current counts, alongside the original nominal constraints. A cap can couple
family sizes even when no other produced/retained practical restriction is
active. Its proxy is an explicitly qualified approximation; optional final
actual-genotype scenarios describe a separate model-based assessment.
\code{\link{st.cross.evaluate}} gives formulas, snapshot inheritance and
metric/uncertainty definitions.

}
\references{
Okoye UJ, Frisch M, Herzog E (2026). Optimization of family sizes in sets
of crosses with a greedy allocation strategy based on automatic
differentiation. Frontiers in Plant Science 17:1727383.
\doi{10.3389/fpls.2026.1727383}.
}
\value{
An object of class \code{st.cross.allocation}, containing:
\describe{
  \item{\code{crosses}}{One row per pooled family: identifiers, produced
  and retained counts, \code{mu}, \code{va}, \code{accuracy},
  \code{intensity}, \code{selected.mean} and pooled produced bounds.}
  \item{\code{contributions.retained}, \code{contributions.produced}}{
  Named nominal parental gene shares for the retained and produced samples.}
  \item{\code{relatedness.retained}, \code{relatedness.produced}}{
  The corresponding \eqn{c'Kc} quantities, or \code{NA} without inherited
  \code{K}.}
  \item{\code{status}}{Search and numerical qualification.}
  \item{\code{objective}}{Predicted retained mean \code{value}, weighted
  \code{response}, fixed \code{unselected.mean}, objective upper
  \code{bound}, response \code{gap}, and \code{bound.type}.}
  \item{\code{constraints}}{Budget, retained total, pooled family count,
  nominal contribution sample, inherited constraints and re-evaluated
  \code{diversity.loss} when an inherited family-slot constraint exists.}
  \item{\code{diagnostics}}{Resource estimates, actual transitions,
  completion state, pooling count and numerical/ancestry qualifications.}
  \item{\code{inputs}}{Original plan and scores, plus \code{input.family}
  mapping original input rows to pooled output families.}
  \item{\code{control}}{Resolved control values.}
}
}
\seealso{\code{\link{st.cross.plan}}, \code{\link{gs.cross.eval.es}},
\code{\link{st.cross.progeny}}, \code{\link{st.options}}}
