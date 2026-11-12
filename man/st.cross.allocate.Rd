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
                  accuracy = 1, control = list())

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
  \item{control}{Named list with unique supported components:
  \describe{
    \item{\code{max.cells}}{Positive integer limit on the estimated DP
    utility, traceback and numeric-row cells; default \code{10000000}.
    This is a DP-workspace estimate, not a total byte-memory limit.}
    \item{\code{max.transitions}}{Positive integer limit on a conservative
    DP transition-count estimate; default \code{100000000}.}
    \item{\code{tolerance}}{Nonnegative response-gap tolerance; default
    \code{1e-10}. A completed search is qualified optimal when its
    enclosed response gap is at most this value times
    \code{max(1, abs(response))}.}
  }}
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
The function uses the same finite-family Burrows selection-intensity
calculation as \code{gs.cross.eval.es(N=r_e, G=g_e)}. It uses zero intensity
when all produced progeny are retained. Accuracy scales the selection
response rather than the family mean.

Repeated reciprocal or identically oriented pair rows describe a single
biological selection family. Their retained counts and produced lower and
upper bounds are added before evaluating selection intensity. Selection
therefore retains the pooled total from the pooled family. It does not
perform independent selection in the original repeated rows.

A two-row integer budget dynamic program enumerates all feasible family
sizes; integer traceback recovers the allocation. The constant weighted
family-mean component is excluded from DP comparisons. No concavity or
continuous-size assumption is required. The algorithm has time complexity
\eqn{O(FBW)} and traceback space \eqn{O(FB)}, where \eqn{F} is the number
of pooled families, \eqn{B} the budget above summed minima, and \eqn{W}
the largest feasible size range. These are fixed-pair allocation bounds,
not claims of a joint optimum over parent pairs and progeny sizes.

When \code{plan} is an \code{st.cross.plan} result, nominal retained
parental contributions must satisfy its inherited contribution bounds,
fixed-contribution requirements and relatedness limit. A retained family
contributes half its nominal gene share to each parent. The inherited
\code{K} is used to report retained and produced group relatedness
separately. Produced contributions are descriptive and are not constrained
by retained-share limits. A family-slot diversity constraint inherited
from the original fixed-pair plan is unchanged by altering produced sizes.
A plain data frame has no inherited \code{K} or planning constraints, so
its relatedness values are \code{NA}.

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
\code{\link{st.options}}}
\examples{
data("v-tropmaize-vcf")
data("v-tropmaize-phe")
st.load.vcf.data(v.tropmaize.vcf, data.set = "allocation.example")
st.load.performance.data(v.tropmaize.phe, data.set = "allocation.example")
st.restrict.marker.data(NoAll.MAX = 2, data.set = "allocation.example")
st.restrict.marker.data(MaMis.MAX = 0.05, data.set = "allocation.example")
st.restrict.marker.data(ExHet.MIN = 0.1, data.set = "allocation.example")
st.restrict.marker.data(InMis.MAX = 0.1, data.set = "allocation.example")
gs.esteff.rr(method = "BLUP", data.set = "allocation.example")
effects <- gs.return.effects(data.set = "allocation.example")
parents <- head(rownames(st.get.A("allocation.example")), 6)
st.restrict.marker.data(ind.list = parents, data.set = "allocation.example")
st.restrict.marker.data(MaMis.MAX = 0, data.set = "allocation.example")
map <- st.get.map(data.set = "allocation.example")
effects <- effects[effects$marker \%in\% c("mu", map$Name), ]
gs.set.effects(effects, data.set = "allocation.example")
gs.cross.eval.mu(data.set = "allocation.example")
gs.cross.eval.va(pop.type = "DH", data.set = "allocation.example")
scores <- gs.cross.info(data.set = "allocation.example")
pairs <- scores[1:3, c("P1Name", "P2Name")]
allocation <- st.cross.allocate(pairs, scores, total.progeny = 60,
                                n.selected = 5, max.progeny = 40)
allocation
allocation$crosses
}
