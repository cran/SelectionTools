\name{st.cross.plan}
\alias{st.cross.plan}
\alias{print.st.cross.plan}
\alias{st.cross.progeny}
\alias{st.cross.frontier}
\alias{print.st.cross.frontier}
\title{Cross Planning, Progeny Budgets and Gain-Diversity Frontiers}
\description{
Optimum cross selection chooses distinct non-selfed parent pairs with, by
default, equal family sizes, contribution
bounds and a group-relatedness limit. Joint planning maximizes additive mean;
an optional cross-score table also supports usefulness-based joint planning.
Fixed-contribution planning minimizes average relatedness within crosses.
The optional \code{progeny} group jointly chooses pairs and produced counts,
with either fixed family retention or an explicit total retained budget.

Cross-planning frontiers report assessed plans, feasibility and numerical
status at a declared sequence of diversity limits. They support both
current-generation and additional-cycle population objectives.
}
\usage{
st.cross.plan(merit, K, n.crosses, theta, n.progeny = 100,
              lower = 0, upper = 0.5, pairs = NULL,
              contributions = NULL, mode = "joint", control = list(),
              scores = NULL, criterion = "mean", accuracy = 1,
              diversity = NULL, constraints = NULL, progeny = NULL,
              genetic.model = NULL)

\method{print}{st.cross.plan}(x, digits = NULL, top = NULL, bottom = NULL,
                            char.per.line = NULL, ...)

st.cross.frontier(merit, K, n.crosses, theta, n.progeny = 100,
                  lower = 0, upper = 0.5, pairs = NULL,
                  contributions = NULL, mode = "joint", control = list(),
                  scores = NULL, criterion = "mean", accuracy = 1,
                  diversity = NULL, constraints = NULL, progeny = NULL,
                  genetic.model = NULL, limits, metric = "theta",
                  warmstart = FALSE)

\method{print}{st.cross.frontier}(x, digits = NULL, top = NULL,
                                bottom = NULL, char.per.line = NULL, ...)
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
  the character labels \code{"1"}, \code{"2"}, and so on. With practical
  \code{constraints}, a full named reference may include additional unused
  parents; named candidate merits identify the planning subset. The full
  reference remains in \code{inputs$K}, and unused reference shares are zero.
  The default role-free matrix/order contract remains unchanged.}
  \item{n.crosses}{Positive supported whole number of distinct families.
  Each unordered parent pair occurs at most once, and self crosses are
  excluded.}
  \item{theta}{Finite upper limit on group relatedness, \eqn{c'Kc}, on the
  original reference scale of \code{K}.}
  \item{n.progeny}{Positive supported whole number of progeny per family.
  All families have this same size on the equal-family path.}
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
  same allowed pair under the default role-free model. With explicit
  reproductive roles in \code{constraints}, column one identifies the seed
  parent and column two the pollen parent; reciprocal directions are
  alternatives for one biological family.}
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
    start, including budgeted degree-realization attempts and local changes;
    default \code{20000}.}
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
    not a continuation at an earlier saved random-number state. The
    \code{progeny} path also accepts an allocation as a count candidate and
    independently rescores it under current inputs.}
  }}
  \item{scores}{Optional data frame with character columns \code{P1Name}
  and \code{P2Name}, and numeric columns \code{mu}, \code{va}, and
  \code{es}. These contain cross means, nonnegative segregation variances,
  and expected superior-progeny scores, as returned by \code{gs.cross.info}
  after the corresponding evaluations. All score values are finite.
  Each allowed pair has exactly one score row. Extra rows are permitted
  but all rows contain known, different parents and unique unordered pairs.
  Row order and parental order within a pair are immaterial. Additional
  columns are permitted. The table is required for usefulness planning.
  With \code{progeny}, only \code{mu} and \code{va} are required and used;
  the selected score is recomputed for actual counts.}
  \item{criterion}{Either \code{"mean"}, \code{"usefulness"}, or
  \code{"future"}.
  Without \code{scores}, mean planning uses the mean of parental merits.
  With \code{scores}, mean and usefulness planning use the supplied cross
  means, and
  usefulness includes the supplied expected within-family selection gain on
  the equal-family path. With \code{progeny}, usefulness uses the requested
  count-dependent selection calculation. Future planning requires an
  explicit \code{progeny} group and genetic model, and targets the actual
  next selected DH population described in \code{\link{st.cross.future}}.}
  \item{accuracy}{Finite within-family selection accuracy between zero
  and one. The usefulness score is
  \code{mu + accuracy * (es - mu)}. The default of one retains the
  supplied \code{es}; zero gives the supplied \code{mu}. A nondefault
  value requires \code{scores}. With \code{progeny}, accuracy multiplies
  the recomputed finite-family selection differential.}
  \item{diversity}{Optional list with exactly \code{G}, \code{max}, and
  \code{sample = "family.slots"}. The sample declaration is mandatory.
  \code{G} is the complete named fixed reference covariance matrix accepted
  by \code{propsd}; it contains every planning parent and may contain
  additional reference parents. \code{max} is a finite strict upper limit
  on the computed proportional loss. \code{NULL} retains the unconstrained
  interface and output structure.}
  \item{genetic.model}{Optional additive DH genetic model requiring an
  explicit \code{progeny} group, as described in
  \code{\link{st.cross.evaluate}}. The model supplies approximate
  selection-aware founder contributions and optional final actual-genotype
  scenarios. Its optional \code{future} group supports actual next-cycle
  assessment and planning as described in \code{\link{st.cross.future}}.
  A supplied frozen snapshot supports use without a live marker data set.
  Planner starts supply counts and pairs, without an implicit genetic
  model or future policy. Scenario execution always requires a current
  request. An active
  \code{selected.founder.theta} requires the original supplied
  \code{K} and couples the actual produced and retained counts.}
  \item{x}{A cross-planning result or cross-planning frontier result.}
  \item{constraints}{Optional practical parent, pool, role and resource
  restrictions described in \code{\link{st.cross.constraints}}. They
  intersect the existing contribution, pair, relatedness and diversity
  restrictions. Retained-count and operation-count restrictions require
  identified counts; the equal-family path does not infer them.}
  \item{progeny}{Optional joint family-size specification described in
  \code{\link{st.cross.progeny}}. \code{NULL} retains the equal-family
  planner. The group supplies a produced budget, fixed common retained count
  or variable-retention total, supported family bounds, selection method
  and allocation controls. It requires \code{mode="joint"} and keyed \code{mu}/\code{va}
  scores. An explicit top-level \code{n.progeny} conflicts with this group.}
  \item{digits, top, bottom, char.per.line}{Optional printing controls.
  Unspecified values follow the corresponding \code{st.options} settings.}
  \item{...}{Additional printing arguments.}

\item{limits}{Nonempty numeric vector of finite declared limits, in the
requested point order. Selected-founder proxy limits are nonnegative.}
\item{metric}{Either \code{"theta"} for original family-slot group
relatedness or \code{"selected.founder.theta"} for the selection-aware
founder proxy. Each point replaces the parameter named by \code{metric};
independent restrictions continue to intersect it. The latter metric
requires a valid genetic model and original relationship reference.
Neither is an actual future-genomic bound.}
\item{warmstart}{Logical scalar. The default \code{FALSE} uses each
point's standalone starting rules and seed. With \code{TRUE}, a checked
previous plan can supply actual pairs and counts as a start; acceptance,
rejection and fallback remain explicit. Held-out results do not choose
warm starts.}
}
\details{
The following equal-family score formulas describe \code{progeny=NULL}.
The optional joint produced-budget model and its count-dependent formulas
are described in \code{\link{st.cross.progeny}}. Relatedness and top-level
parental bounds retain their explicit family-slot basis on that path.

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
parental use counts. A successful greedy construction or feasible warm start
proceeds directly to local search. If construction fails before a feasible
incumbent is available, randomized degree realization can use the remaining
proposal budget for that start. Each attempt counts against \code{control$maxit};
these attempts provide no optimality certificate. Failure to find a feasible
plan within a search budget is distinguished from proven infeasibility.

Mean and usefulness describe the declared additive model. Parental
contributions remain nominal shares before within-family selection.
Realized ancestry after selecting descendants
is not determined by these contributions alone. An equal-family plan can be
executed with \code{st.cross.plan.execute} when the selected founders satisfy
that function's simulation requirements. Varying actual family counts are
available to the individual simulation commands.

Inputs, the caller's R and simulation random-number states, and thread
settings are preserved by the planner. The local search seed and final state
are returned with the result. The numeric input precision is retained.
The optional genetic model and its inherited selected-founder cap use actual
current counts, alongside the original nominal constraints. A cap can couple
family sizes even when no other produced/retained practical restriction is
active. Its proxy is an explicitly qualified approximation; optional final
actual-genotype scenarios describe a separate model-based assessment.
\code{\link{st.cross.evaluate}} gives formulas, snapshot inheritance and
metric/uncertainty definitions.

}
\section{Joint progeny budgets}{
The named group has the following supported components:
\describe{
  \item{\code{total.progeny}}{Required positive supported whole-number
  produced budget, at most the C integer maximum.}
  \item{\code{n.selected}}{Positive whole-number scalar retained count per
  active family for the fixed-retention contract. Total retention is
  \code{n.crosses*n.selected}. It is mutually exclusive with the variable
  retained-budget components below.}
  \item{\code{total.selected}}{Required positive supported whole-number
  retained budget for variable retention, at most \code{total.progeny}.
  Actual family counts sum to this budget.}
  \item{\code{min.progeny}}{Common positive whole-number lower bound,
  defaulting to \code{n.selected} under fixed retention and to one under
  variable retention. Every actual produced count is at least its actual
  retained count.}
  \item{\code{max.progeny}}{Common upper bound, defaulting to the total
  produced budget, no smaller than the lower bound.}
  \item{\code{min.selected}, \code{max.selected}}{Common retained-count
  bounds for variable retention. The lower bound defaults to one and can
  explicitly be zero; the upper defaults to the smaller of the resolved
  common \code{max.progeny} and \code{total.selected}. Both are nonnegative
  whole numbers, with the upper no smaller than the lower. Actual retained
  counts are additionally bounded by actual production.}
  \item{\code{bounds}}{Optional variable-retention data frame with character
  keys \code{P1Name}, \code{P2Name} and one or more of \code{min.progeny},
  \code{max.progeny}, \code{min.selected}, \code{max.selected}. Each
  unordered biological family occurs at most once and belongs to the
  eligible graph. Common scalar defaults resolve first; provided columns
  override them for the keyed family. Missing rows and absent columns inherit
  those defaults. Missing cells, unknown columns and duplicate keys produce
  errors. Reciprocal directions share one biological count domain.
  These bounds require \code{total.selected}.}
  \item{\code{selection.method}}{\code{"burrows"} by default, or
  \code{"order.statistics"}; the latter uses the qualified finite-normal
  order-statistic calculation of \code{st.cross.allocate}.}
  \item{\code{control}}{Named allocation controls \code{max.cells}
  (default \code{10000000}), \code{max.transitions} (\code{100000000}),
  \code{tolerance} (\code{1e-10}) and \code{algorithm} (\code{"auto"}).
  Workspace includes shared caches and live inner-solver work. The transition
  allowance covers the whole joint call, rather than restarting for every
  proposed pair set. The group accepts \code{"auto"} or \code{"dp"};
  \code{"greedy"} additionally requires \code{"order.statistics"}.
  Coupled practical count restrictions require \code{"auto"}.
  Variable retention supports \code{"auto"} and the compatible separable
  \code{"dp"} request; it does not support \code{"greedy"}.
  Unsupported combinations produce errors.}
}

The top-level planner \code{control} still supplies its outer search method,
seed, starts, proposals and node cap. The group is supported in joint mode.
An explicitly supplied top-level \code{n.progeny} together with this group is
an error; its default applies only to the equal-family path.
Candidate scores require finite \code{mu} and nonnegative finite \code{va}
for each allowed unordered pair. Stored \code{es} does not determine the
count-dependent selection response.

With \eqn{C} active families, produced counts \eqn{n_e}, common retained
count \eqn{r}, accuracy \eqn{\rho} and intensity \eqn{a(n_e,r)},
\deqn{\sum_e n_e=P,\qquad
U_e=\mu_e+\rho\sqrt{v_e}\,a(n_e,r),\qquad
J_{G1}=C^{-1}\sum_{e\ \mathrm{active}}U_e.}
These formulas describe fixed common retention. With a total retained budget
\eqn{R}, actual retained counts \eqn{r_e} instead satisfy
\deqn{\sum_e r_e=R,\quad 0\leq r_e\leq n_e,\qquad
J_{G1}=\sum_{e:r_e>0}\frac{r_e}{R}
 \{\mu_e+\rho\sqrt{v_e}\,a(n_e,r_e)\}.}
The corresponding mean criterion is \eqn{\sum_e(r_e/R)\mu_e}.
Thus retention changes the weights of cross means as well as their selection
intensities. Selection is still within families; a pooled top-\eqn{R} rule
over all produced descendants is a different model.
An inactive family produces and retains zero progeny. At most one eligible
direction is selected for each unordered biological family.
The usefulness criterion maximizes \eqn{J_{G1}}; the mean criterion maximizes
the unselected family-mean baseline on the same retained weights.

Every active family produces at least one descendant and consumes one of
the requested family slots. With an explicit zero retained lower bound, an
active family can retain none. It still contributes to produced and family
metrics, while its retained weight is zero. Its selected mean and intensity
are unavailable; no intensity at zero retention is evaluated. Totals skip
zero-weight rows explicitly. Full retention gives exactly zero selection
response. The number of families contributing retained descendants can be
smaller than the requested produced-family count.

Top-level parental contribution bounds, \code{theta} and family-slot
diversity retain their family-slot basis. Practical restrictions additionally
use their explicitly declared produced or retained bases. Original named
\code{K}/\code{G} references and scales remain fixed. Unknown operation
counts are not inferred from produced counts.

Exact small search enumerates biological pair/direction choices and solves
the corresponding supported fixed-pair allocation, with valid bounds for
the requested objective. Larger searches use bounded starts and exchanges.
An exact inner allocation does not establish an outer joint optimum.
An unfinished inner search leaves an unresolved branch unless an independent
bound excludes it. Exhausting search or repair budgets does not prove
infeasibility. A feasible supplied incumbent remains available when limits
stop further improvement.

Variable-retention exact small solving uses both production and retention
budgets. Larger problems retain checked incumbents with bounded count and
family changes. The decreasing-produced-increment argument for fixed
retention does not establish a joint optimum when retention changes.
Unweighted family-average pruning does not bound the retained-weighted
objective. A valid weighted relaxation or an unknown bound is reported.
Original means are centered, with baseline and response compared separately
in mean units, to preserve small responses under large common merit offsets.

The count-only intensity and discrete marginal cache is shared across the
joint call. A change in parental merit or variance rescales the cached
kernel rather than changing its mathematical definition. Numerical normal
quadrature error estimates remain estimates rather than rigorous intervals.
Complete integer search and biological prediction accuracy are distinct
qualifications.

Produced-budget results have class \code{st.cross.plan}, with actual produced/retained
counts, recomputed selected means and explicit method/basis provenance.
\code{st.cross.evaluate} evaluates the same supplied population.
A subsequent \code{st.cross.allocate} intentionally supplies new production
bounds/budget for its fixed pairs, while retaining the source plan's
operational and parental restrictions. It does not reimpose the source
equal-size design used as a baseline.

Individual simulation commands accept the actual counts of varying families.
The separate \code{st.cross.plan.execute} interface retains its equal-count
execution contract.

A saved feasible plan is a warm start through \code{control$start}; current
inputs independently rescore and check it. A warm restart is a new search
and does not restore an unrecorded exact-search frontier.
Variable-retention warm starts require actual \code{n.selected} counts;
divisibility of the retained budget by the family count does not identify
those counts. A fixed-retention plan can supply a start when its actual
counts satisfy both new budgets and all current family bounds.
}
\section{Gain-diversity frontiers}{
A common immutable genetic core, simulation banks and reusable caches
support all frontier points. Genetic cell/work limits are cumulative
across the whole call, including final assessments; forward simulation does
not receive a new budget at each point.

Completed points reproduce standalone schedules with the recorded inputs,
seed, start and bank when sufficient matching budgets are available.
A point stopped by the shared budget is explicitly incomplete and need
not match an uncapped standalone search. Infeasible, unavailable and
unresolved outcomes remain in the requested sequence.

The frontier is a collection of assessed solutions, not a certificate
between points. Heuristic or finite-bank points need not be monotone.
Current family-slot quantities, approximate selected-founder proxies and
actual current/future genomic metrics retain separate labels. Future
assessment and paired uncertainty follow \code{\link{st.cross.future}}.
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
The table has zero rows when no feasible plan is returned. With
\code{progeny}, the table contains actual \code{n.selected}, selection
\code{accuracy}, \code{intensity}, \code{selected.mean} and \code{response};
stored \code{es} does not supply its count-dependent utility.}
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
includes completed enumeration with unresolved boundary feasibility.
The finite-normal joint produced-budget mode separately qualifies estimated
quadrature bounds; numerical optimality does not constitute a rigorous
interval certificate or validation of the biological normal approximation.
Its additional statuses and cumulative-work diagnostics are described in
\code{\link{st.cross.progeny}}.}
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
incumbent was used. The proposal count includes budgeted degree-realization
attempts as well as local-change proposals.
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

For the \code{progeny} argument, the returned plan
identifies the solved count-dependent model, feasible incumbent, bound/gap
and numerical/search qualification. Actual \code{crosses$n.progeny} and
\code{crosses$n.selected} define its resources.

\code{optimal} qualifies the Burrows objective when its gap is closed;
\code{numerically.optimal} additionally identifies the estimated numerical
qualification of normal order statistics. A checked feasible result with
an open or unknown gap has \code{best.found} status. Numerical boundary
uncertainty can produce \code{numerically.unresolved}; absence of a feasible
incumbent under search or allocation limits gives \code{search.limit} or
\code{resource.limit}. Only a qualified impossibility result is
\code{infeasible}. Shared cache, work, workspace and timing diagnostics
are reported separately from the mathematical objective.

An object of class \code{st.cross.frontier} containing the requested
limits and metric, point plans and statuses, summaries, starting provenance
and shared resource diagnostics. Printing returns the input invisibly.
}
\seealso{
\code{propsd}, \code{st.ocs}, \code{st.cross.plan.execute}, \code{gs.cross.info},
\code{gs.cross.eval.mu}, \code{gs.cross.eval.va}, \code{gs.cross.eval.es},
\code{gs.estimate.gv}, \code{st.get.A}, \code{st.options},
\code{\link{st.cross.progeny}}, \code{\link{st.cross.evaluate}},
\code{\link{st.cross.frontier}}, \code{\link{st.cross.future}}
}
