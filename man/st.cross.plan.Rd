\name{st.cross.plan}
\alias{st.cross.plan}
\alias{print.st.cross.plan}
\title{Optimum Cross Selection with Integer Families}
\description{
Optimum cross selection chooses distinct non-selfed parent pairs with equal
family sizes, contribution
bounds and a group-relatedness limit. Joint planning maximizes additive mean;
an optional cross-score table also supports usefulness-based joint planning.
Fixed-contribution planning minimizes average relatedness within crosses.
}
\usage{
st.cross.plan(merit, K, n.crosses, theta, n.progeny = 100,
              lower = 0, upper = 0.5, pairs = NULL,
              contributions = NULL, mode = "joint", control = list(),
              scores = NULL, criterion = "mean", accuracy = 1,
              diversity = NULL)

\method{print}{st.cross.plan}(x, digits = NULL, top = NULL, bottom = NULL,
                            char.per.line = NULL, ...)
}
\arguments{
  \item{merit}{Finite numeric vector of additive parental merits, with at
  least two parents. Its order agrees with both dimensions of \code{K}.
  Names, when supplied, are unique, complete parental IDs.}
  \item{K}{Finite symmetric positive-semidefinite relationship matrix in
  parental order. Singular matrices are supported. For an additive
  relationship matrix \code{A}, the kinship convention is \code{K = A/2}.
  Row and column names are either both absent or both supplied and equal.
  Named \code{merit} entries agree with these names in the same order.
  Parental IDs follow \code{merit} names, otherwise \code{K} names, otherwise
  the character labels \code{"1"}, \code{"2"}, and so on.}
  \item{n.crosses}{Positive supported whole number of distinct families.
  Each unordered parent pair occurs at most once, and self crosses are
  excluded.}
  \item{theta}{Finite upper limit on group relatedness, \eqn{c'Kc}, on the
  original reference scale of \code{K}.}
  \item{n.progeny}{Positive supported whole number of progeny per family.
  All families have this same size.}
  \item{lower, upper}{Lower and upper bounds on parental gene contributions.
  Each is an unnamed numeric scalar or a numeric vector matching
  \code{merit}. Full-length names, when supplied, match the parental order.
  Bounds satisfy \code{0 <= lower <= upper <= 1}; a non-selfed plan can use
  no parent in more than half of its total parental gene contribution.}
  \item{pairs}{Allowed unordered pairs, supplied as a two-column character
  matrix or a data frame with exactly two character columns containing
  parental IDs. Column order specifies the two parents; column names are not
  used. Unknown IDs, self pairs and repeated unordered pairs produce errors.
  \code{NULL} allows every non-selfed pair. Reciprocal rows describe the
  same allowed pair.}
  \item{contributions}{Required only for \code{mode = "fixed"}. A finite
  numeric vector in parental order specifying nonnegative gene contributions
  that sum to one. Full-length names, when present, match the parental IDs.
  Each contribution must be realizable as an integer number of parental uses
  divided by \code{2 * n.crosses}, and the resulting uses must satisfy the
  contribution bounds and allowed-pair restrictions.}
  \item{mode}{Either \code{"joint"}, which jointly chooses parental usage
  and pairing under the selected criterion, or \code{"fixed"}, which
  preserves supplied contributions and minimizes average within-cross
  relatedness. Fixed mode requires \code{scores = NULL} and
  \code{criterion = "mean"}.}
  \item{control}{Named list of search controls. Unknown and repeated control
  names produce errors. The supported components are:
  \describe{
    \item{\code{method}}{\code{"auto"} by default, or \code{"exact"}
    or \code{"heuristic"}. Automatic selection uses the exact search when
    its estimated work fits the node budget, and the heuristic otherwise.}
    \item{\code{max.nodes}}{Positive whole-number exact-search node budget;
    default \code{1000000}.}
    \item{\code{maxit}}{Positive whole-number proposal budget per heuristic
    start; default \code{20000}.}
    \item{\code{starts}}{Positive whole number of heuristic starts;
    default \code{16}.}
    \item{\code{seed}}{Nonnegative supported whole-number seed for the
    search's local random-number stream; default \code{1}.}
    \item{\code{tolerance}}{Finite number strictly between zero and one;
    default \code{1e-8}. Controls numerical comparisons. The relatedness
    allowance uses the active quadratic terms and requested limit, as
    described below. Reported objective and relatedness values use the
    original input units.}
    \item{\code{start}}{Optional earlier \code{st.cross.plan} result or
    its \code{crosses} data frame, containing a feasible incumbent with the
    requested family count and progeny count, whose pairs satisfy the current
    allowed-pair restrictions and constraints. Its actual family rows are
    validated before use. This starts a new search from an incumbent; it is
    not a continuation at an earlier saved random-number state.}
  }}
  \item{scores}{Optional data frame with character columns \code{P1Name}
  and \code{P2Name}, and numeric columns \code{mu}, \code{va}, and
  \code{es}. These contain cross means, nonnegative segregation variances,
  and expected superior-progeny scores, as returned by \code{gs.cross.info}
  after the corresponding evaluations. All score values are finite.
  Each allowed pair has exactly one score row. Extra rows are permitted
  but all rows contain known, different parents and unique unordered pairs.
  Row order and parental order within a pair are immaterial. Additional
  columns are permitted. The table is required for usefulness planning.}
  \item{criterion}{Either \code{"mean"} or \code{"usefulness"}.
  Without \code{scores}, mean planning uses the mean of parental merits.
  With \code{scores}, both criteria use the supplied cross means, and
  usefulness includes the supplied expected within-family selection gain.}
  \item{accuracy}{Finite within-family selection accuracy between zero
  and one. The usefulness score is
  \code{mu + accuracy * (es - mu)}. The default of one retains the
  supplied \code{es}; zero gives the supplied \code{mu}. A nondefault
  value requires \code{scores}.}
  \item{diversity}{Optional list with exactly \code{G}, \code{max}, and
  \code{sample = "family.slots"}. The sample declaration is mandatory.
  \code{G} is the complete named fixed reference covariance matrix accepted
  by \code{propsd}; it contains every planning parent and may contain
  additional reference parents. \code{max} is a finite strict upper limit
  on the computed proportional loss. \code{NULL} retains the unconstrained
  interface and output structure.}
  \item{x}{A cross-planning result.}
  \item{digits, top, bottom, char.per.line}{Optional printing controls.
  Unspecified values follow the corresponding \code{st.options} settings.}
  \item{...}{Additional printing arguments.}
}
\details{
For \eqn{m} equal-sized families, a parent used in \eqn{d_i} families has
nominal gene contribution \eqn{c_i = d_i/(2m)}. Each selected pair contributes
one use to each of its two parents. The contributions sum to one. Equal
positive family sizes make family shares and produced-progeny shares equal;
changing the common family size does not change the selected-pair objective.

In joint mode the objective is
\deqn{
\frac{1}{m}\sum_{(i,j)\in S}\frac{g_i+g_j}{2} = g'c,
}{
mean over selected pairs of (g[i]+g[j])/2 = sum(g*c),
}
where \eqn{S} contains the selected pairs and \eqn{g} is \code{merit}.
The planner chooses both usage and pairing under \eqn{c'Kc\leq\theta} and
the contribution bounds. Integer usage bounds follow from the supplied
fractional bounds, without allowing extra uses beyond those bounds.

When \code{scores} is supplied, joint mean planning instead maximizes the
average supplied \code{mu} over selected pairs. Usefulness planning
maximizes the average
\deqn{U_{ij} = \mu_{ij} + r(es_{ij}-\mu_{ij}),}{
U[i,j] = mu[i,j] + accuracy * (es[i,j] - mu[i,j]),
}
where \eqn{r} is \code{accuracy}. Both criteria retain the same constraints
and report both mean and usefulness for their selected pairs. The supplied
cross means need not equal means reconstructed from exported parental
merits. \code{parent.gain} separately reports \eqn{g'c} in this case.

The upstream functions \code{gs.cross.eval.mu}, \code{gs.cross.eval.va},
\code{gs.cross.eval.es}, and \code{gs.cross.info} supply one common score
table for comparisons. \code{gs.cross.eval.es(N = s, G = n)} incorporates
the existing finite-family approximation of Burrows (1972); its \code{N}
and \code{G} are retained and produced individuals per family.
The planner does not recalculate this score when \code{n.progeny} changes.
The family size, retained count, progeny type, marker-effect model and
selection convention represented by the table therefore form part of the
declared analysis. In particular, usefulness is a fixed cross score in this
interface, rather than an allocation of varying family sizes.

The table requires \code{es >= mu}, apart from arithmetic-precision noise;
\code{control$tolerance} does not relax score validation. Supplied score
values are retained without rounding or replacement. An accuracy of zero,
or an upstream score with \code{es = mu}, reduces usefulness planning to
mean planning on the same table. Fixed cross scores leave the objective
linear in pair selection; they do not change the interpretation of the
group-relatedness constraint.

In fixed mode the objective is
\deqn{
\frac{1}{m}\sum_{(i,j)\in S}K_{ij}.
}{
mean over selected pairs of K[i,j].
}
Both \eqn{g'c} and \eqn{c'Kc} are fixed by the supplied contributions.
Re-pairing those contributions can reduce average within-cross relatedness,
but cannot change the overall additive mean. Arbitrary continuous
contributions from \code{st.ocs} need not correspond to an integer plan.
Machine-arithmetic noise near an integer use count is tolerated; other
nonintegral requests are reported as infeasible rather than silently rounded.
Increasing \code{control$tolerance} does not change this integer-realizability
requirement. Integer
uses also need a realization with distinct allowed non-selfed pairs.

The matrix and its allele-frequency reference remain fixed throughout the
search. A centered genomic matrix describes a relatedness index relative to
that reference. Dividing an additive relationship matrix by two does not
itself provide an absolute probability of identity by descent. The
\code{relatedness} limit and the average \code{within.cross.relatedness}
are separate quantities.

Relatedness is recomputed from the original \code{K} and integer parental
uses. Its numerical feasibility allowance is
\deqn{
\mathrm{tolerance}\,\max\left(|\theta|,
  \sum_{i,j}|K_{ij}|c_i c_j\right).
}{
tolerance * max(abs(theta), sum(abs(K[i,j])*c[i]*c[j])).
}
Matrix entries for unused parents do not create an extra allowance. The
entire computed roundoff interval must satisfy the bound including this
allowance; an interval that crosses its boundary remains unresolved. An exact
zero limit additionally requires an arithmetic-precision check of the
original matrix rows applied to \eqn{c}. This zero check is independent of
\code{control$tolerance}. Numerically unresolved boundary points do not
establish infeasibility or optimality; a result may retain
\code{"search.limit"} or \code{"best.found"} for this reason.

With \code{diversity} enabled, each chosen family contributes one occurrence
of each parent to the selected multiset used by \code{propsd}. Thus parent
\eqn{i} has multiplicity \eqn{d_i} and there are \eqn{2m} occurrences.
This family-slot convention is explicit and does not represent produced
progeny counts, retained descendant counts or continuous gene shares. The full
reference matrix is retained, including unused reference parents. The
constraint is the strict computed comparison \code{loss < max}; the search
comparison tolerance does not relax it. This constrains a Gamma approximation,
not a formal interval enclosure of biological diversity loss. The statistic
can be negative and is not clipped. It remains separate from both relatedness
quantities and from the mean or usefulness objective.

The exact search checks diversity at each complete candidate; its unconstrained
objective bounds remain valid but may be loose. The heuristic constructions
can initially be infeasible. In joint mode, repair uses both relatedness and
diversity violations, and a move from a feasible state must preserve both
computed constraints. Only independently rechecked feasible incumbents are
retained. In fixed mode, multiplicities and therefore PropSD are fixed, so a
violated diversity limit establishes infeasibility before rewiring. Failure to
find a feasible joint plan heuristically remains a search limit.

The exact method uses a bounded exhaustive search with pruning. An optimal
status requires a completed search with a feasible optimum certified under
the stated floating-point comparisons. A search stopped by the node budget
does not receive that status. The
heuristic uses multiple constructions and local changes while retaining only
independently rechecked feasible plans. In fixed mode, local rewiring preserves
parental use counts. Failure to find a feasible plan within a search budget
is distinguished from proven infeasibility.

Mean and usefulness describe the declared additive model. Parental
contributions remain nominal shares before within-family selection.
Realized ancestry after selecting descendants
is not determined by these contributions alone. The plan can be executed
with \code{st.cross.plan.execute} when the selected founders satisfy that
function's simulation requirements.

Inputs, the caller's R and simulation random-number states, and thread
settings are preserved by the planner. The local search seed and final state
are returned with the result. The numeric input precision is retained.
}
\value{
A list of class \code{st.cross.plan}. Its printed representation summarizes
the objective, feasibility, families and parental usage under the current
printing options. The principal components are:
\item{crosses}{Data frame with columns \code{family}, \code{P1Name},
\code{P2Name}, \code{n.progeny}, \code{mu}, and \code{kinship}.
Without a score table, \code{mu} is the mean of the two supplied parental
merits. With a score table it is the supplied cross mean, and the table
also contains \code{va}, \code{es}, and \code{usefulness}.
\code{kinship} is the pair's entry in the supplied \code{K}.
The table has zero rows when no feasible plan is returned.}
\item{contributions}{Named numeric vector of achieved parental gene shares,
or \code{NULL} when no feasible plan is returned.}
\item{parents}{Data frame in parental order, containing \code{parent},
\code{uses}, \code{contribution}, \code{lower}, \code{upper}, and
\code{merit}, including unused parents.}
\item{gain}{Achieved mean merit: \eqn{g'c} without a score table, or the
average supplied cross mean with a score table.}
\item{usefulness, parent.gain, criterion, accuracy}{Additional components
when a score table is supplied: average achieved usefulness, parental
merit \eqn{g'c}, the requested criterion and within-family accuracy.}
\item{diversity}{Present only when enabled. Contains the requested \code{max}
and \code{sample}, achieved \code{loss}, a \code{feasible} indicator, and
\code{evaluation}, the complete \code{propsd} result. When no feasible plan
is returned, loss is \code{NA}, feasibility is false and evaluation is
\code{NULL}. The request is also retained in \code{inputs$diversity}, with
the threshold and sample declaration in \code{constraints$diversity}.}
\item{relatedness}{Achieved group relatedness, \eqn{c'Kc}.}
\item{within.cross.relatedness}{Mean \eqn{K_{ij}} over the selected pairs.}
\item{mode}{The requested planning mode.}
\item{status}{\code{"optimal"} for a completed exact search with a certified
feasible optimum; \code{"best.found"} for a feasible incumbent without a completed
optimality proof; \code{"infeasible"} for a proved infeasible request;
or \code{"search.limit"} when search or numerical-certification limits
leave no certified feasible plan and infeasibility was not proved. This
includes completed enumeration with unresolved boundary feasibility.}
\item{objective}{Criterion, optimization direction, achieved value,
relaxation bound, remaining gap, and bound description. A joint-mode bound is
an upper bound on the selected mean or usefulness criterion; a fixed-mode
bound is a lower bound on average within-cross relatedness.}
\item{constraints}{Requested family counts, progeny budget, relatedness
limit, contribution bounds, effective integer use bounds, requested fixed
contributions, and the number of allowed pairs.}
\item{diagnostics}{Constraint residuals, infeasibility or search-limit reason,
node and proposal counts, accepted moves, starts, completed-search indicator,
method used, relaxation bounds, estimated work, and whether a starting
incumbent was used.
Boundary diagnostics include \code{relatedness.roundoff}, the roundoff
half-width in the reported relatedness units, including numeric-conversion
uncertainty; \code{zero.row.residual}, the largest row-relative residual
for a zero limit; and \code{uncertain.evaluations}, the number of candidate
evaluations with unresolved boundary feasibility. The reason
\code{"relatedness.numerical.limit"} identifies unresolved relatedness
certification.}
\item{control}{Search controls used for the request.}
\item{search}{Local seed, final random-number state, generator label and
whether a supplied starting plan was used.}
\item{inputs}{The supplied merit vector, relationship matrix and allowed
pairs, plus the optional score table, retained for reconstruction and review.}
}
\seealso{
\code{propsd}, \code{st.ocs}, \code{st.cross.plan.execute}, \code{gs.cross.info},
\code{gs.cross.eval.mu}, \code{gs.cross.eval.va}, \code{gs.cross.eval.es},
\code{gs.estimate.gv}, \code{st.get.A}, \code{st.options}
}
\examples{
\dontrun{
data("v-tropmaize-vcf")
data("v-tropmaize-phe")
st.load.vcf.data(v.tropmaize.vcf)
st.load.performance.data(v.tropmaize.phe)
st.restrict.marker.data(NoAll.MAX = 2)
st.restrict.marker.data(MaMis.MAX = 0.05)
st.restrict.marker.data(ExHet.MIN = 0.1)
st.restrict.marker.data(InMis.MAX = 0.1)
gs.esteff.rr(method = "BLUP")
values <- gs.estimate.gv()
K <- st.get.A()/2
merit <- setNames(values$yhat, as.character(values$i))
merit <- merit[rownames(K)]
parents <- head(names(merit), 6)
merit <- merit[parents]
K <- K[parents, parents, drop = FALSE]
joint <- st.cross.plan(merit, K, n.crosses = 3, theta = max(K),
                       control = list(method = "exact"))
joint
chosen.crosses <- joint$crosses
chosen.crosses
fixed <- st.cross.plan(merit, K, n.crosses = 3, theta = max(K),
                       contributions = joint$contributions, mode = "fixed",
                       control = list(method = "exact", start = joint))
fixed
parent.usage <- fixed$parents
parent.usage
}
}
