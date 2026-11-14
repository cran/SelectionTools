\name{st.cross.evaluate}
\alias{st.cross.evaluate}
\alias{print.st.cross.population}
\alias{st.cross.future}
\encoding{UTF-8}
\title{Evaluation of a Produced and Selected Progeny Population}
\description{
Evaluates the actual integer family sizes in a cross plan without changing
them. Population summaries distinguish produced progeny from retained
progeny and nominal parental contributions on three separate bases.

An optional additive DH genetic model supplies selection-aware founder
proxies and actual-genotype assessments. Its \code{future} group describes
one additional crossing and selection cycle.
}
\usage{
st.cross.evaluate(plan, scores = NULL, n.selected = NULL,
                  accuracy = NULL, K = NULL, selection.method = NULL,
                  constraints = NULL, genetic.model = NULL)

\method{print}{st.cross.population}(x, digits = NULL, top = NULL,
                                  bottom = NULL, char.per.line = NULL, ...)
}
\arguments{
  \item{plan}{An \code{st.cross.plan}, \code{st.cross.allocation}, or data
  frame with character \code{P1Name}, \code{P2Name} and integer-valued
  \code{n.progeny} columns. Unordered repeated families are pooled before
  selection response is evaluated. Self crosses are excluded.}
  \item{scores}{Optional keyed data frame with \code{P1Name},
  \code{P2Name}, finite numeric \code{mu} and nonnegative finite
  \code{va}. Otherwise the plan's stored score columns are used.
  Conflicting scores for repeated families are errors.}
  \item{n.selected}{Nonnegative retained counts, supplied explicitly as
  a scalar or one value per input row, or taken from the stored column.
  Each count is no larger than its produced count.}
  \item{accuracy}{Prediction accuracy in \eqn{[0,1]}, supplied as a scalar
  or one value per input row. A stored column is used when available, followed by a source cross plan's
  stored scalar accuracy; otherwise the value is one. Repeated families have identical accuracy.}
  \item{K}{Optional full named symmetric positive-semidefinite parental
  relationship reference, or the inherited reference. Unused named
  parents have zero contributions. Its original scale is preserved.
  An explicit reference conflicting with inherited constraints is an error.}
  \item{selection.method}{\code{"burrows"} or \code{"order.statistics"}.
  \code{NULL} uses the source allocation's method when stored and the
  source plan's explicit \code{progeny$selection.method} for that model,
  with the Burrows approximation otherwise. The latter choice preserves the
  selection calculation of \code{gs.cross.eval.es}; normal order statistics
  use the finite-family numerical evaluator of \code{st.cross.allocate}.}
  \item{constraints}{Optional practical restrictions described in
  \code{\link{st.cross.constraints}}, intersecting supported inherited
  restrictions. Evaluation checks actual counts and preserves the source
  optimizer's status separately.}
  \item{genetic.model}{Optional additive DH model. A fresh named list has
  \code{data.set}, \code{population="DH"}, \code{t=0},
  \code{map.function="Haldane"},
  \code{index.model="trait.independent.noise"},
  \code{marker.weights=NULL}, \code{selected.founder.theta=NULL},
  \code{scenarios=NULL}, \code{future=NULL} and
  \code{control=list(max.cells=1e7,max.work=1e8)}.
  Only \code{data.set} is required for a fresh model; the other values are
  defaults. Named nonnegative marker weights cover every model marker and
  are normalized once. The model requires complete homozygous biallelic
  candidate parents, complete reference genotypes, fitted additive effects,
  dosage coding 0/1/2 and a complete map in cM. The optional finite nonnegative
  \code{selected.founder.theta} caps the approximate selected-founder
  quadratic on the original supplied or inherited \code{K}.
  \code{scenarios=list(replicates=S,seed=0L)} requests final actual-genotype
  Monte Carlo evaluation, with positive whole \code{S} and a nonnegative
  supported whole-number seed. Genetic work and memory controls are
  cumulative within the call. An explicit \code{future} group requests
  actual next-cycle DH assessment described in \code{\link{st.cross.future}}.
  \code{NULL} inherits an actual source object's
  frozen core and cap, without executing its old scenario request. A list
  without \code{data.set} can request new scenarios, future assessment,
  controls or a tighter
  inherited cap. Conflicting cores and unknown or repeated fields are errors.
  The details below define the model and available summaries.}
  \item{x}{An \code{st.cross.population} result.}
  \item{digits, top, bottom, char.per.line}{Print controls. \code{NULL}
  uses the corresponding \code{st.options()} setting.}
  \item{...}{Reserved for S3 printing.}
}
\details{
For family \eqn{e}, produced and retained counts are \eqn{n_e} and
\eqn{r_e}, with totals \eqn{P} and \eqn{R}. The produced population has
\deqn{M_0=\sum_e(n_e/P)\mu_e,\qquad
V_0=\sum_e(n_e/P)v_e+\sum_e(n_e/P)(\mu_e-M_0)^2.}
The first variance term is within families; the second is between families.
They describe marginal individual variance, not the expected empirical
sample variance.

The expected selected-family mean is
\deqn{U_e=\mu_e+\rho_e\sqrt{v_e}\,a(n_e,r_e),}
where \eqn{a} is the selected finite-family intensity. The selected
population mean is \eqn{J_1=\sum_e(r_e/R)U_e}; its response is retained
separately from the retained-weight baseline \eqn{\sum_e(r_e/R)\mu_e}.
The between-family term \eqn{\sum_e(r_e/R)(U_e-J_1)^2} is reported.
Without an optional genetic model, the selected-population summary contains
the expected mean and its between-family component. It does not provide
selected within-family or total variance, realized genomic ancestry or
selected genomic diversity. Additional-cycle assessment requires the
explicit \code{future} group described below.

Nominal contributions count one half from each parent, weighted either by
active unique families (\code{family.slots}), produced counts
(\code{produced}), or retained counts (\code{retained}). They are not
selection-adjusted genomic ancestry. Relatedness on each available basis
is \eqn{c'Kc}. Produced counts do not identify repeated crossing operations.

Zero-produced families are inactive. A produced family retaining no progeny
has unavailable intensity/selected mean and zero retained weight. Summaries
whose denominators are zero remain unavailable with a reason. Counts and
pooled counts lie between zero and the supported C integer maximum; totals
are checked for exact representation.

Inherited supported budgets, bounds and diversity restrictions are checked
again against actual counts. A structurally valid violation is returned as
infeasible, with residuals; an optimizer's saved status is not proof of
feasibility. Unsupported inherited restrictions remain unqualified.
The evaluator neither reallocates progeny nor certifies optimality.
Numerical order-statistic error estimates are not rigorous interval bounds.

A direct \code{progeny} source plan is checked against its declared active
family count, produced budget and count bounds. Fixed retention supplies a
common retained count; variable retention supplies a total retained budget
and biological-keyed produced/retained family bounds. Zero-retention
produced families remain active on the family and produced bases.
An allocation descended from that plan instead has its own direct production
and retention domain. Original source parental, pairing and practical
restrictions remain inherited. Explicit source contribution and relatedness
bases are preserved; actual produced and retained summaries remain separate.
With \code{genetic.model}, a frozen versioned snapshot stores the candidate
marker genotypes, fitted effects, map, weights and fixed genomic reference.
Subsequent source evaluation does not require the live dataset. Preparation
reads the marker data without changing its datasets, caches or simulation
state. A conflicting explicit \code{K} remains an error; the inherited
original reference and its scale remain fixed. A warm-start plan supplies
pair/count candidates, without implicitly supplying a new planner's model.

Here \eqn{x_{pj}} is parent \eqn{p}'s effect-allele presence and \eqn{u_j}
is the stored dosage coefficient. Parental value is
\eqn{u_0+2\sum_j u_jx_{pj}}. For an immediate DH family,
\eqn{d_j=u_j(x_{1j}-x_{2j})} and \eqn{G=\mu+d'S}, where parental origin
\eqn{S_j} is \eqn{+1} or \eqn{-1}. The origin covariance is
\eqn{Q_{jk}=\exp(-2|m_j-m_k|)} on the same chromosome, with positions
\eqn{m} in Morgans; it is zero between chromosomes. Thus
\eqn{v=d'Qd}, weighted P1 origin is \eqn{C_1=(1+w'S)/2}, and
\eqn{h=\operatorname{Cov}(C_1,G)=w'Qd/2}. Positive linear-time Haldane
recurrences prepare these family moments without a dense marker matrix.

For positive retained count, the selected-origin approximation is
\deqn{c_1=\frac12+\rho\,a(n,r)h/\sqrt{v},\qquad c_2=1-c_1.}
This adapts Allier et al. (2019), G3 equations (4)--(6), and Frontiers in
Genetics equations (5)--(7), using the requested finite-family intensity.
It assumes independent noise in the declared selection index; accuracy alone
does not identify an arbitrary index's ancestry response. Full retention,
zero accuracy and zero trait variance give nominal one-half contributions.
Zero-retention families have no conditional selected contribution and have
zero aggregate weight. The retained-weight aggregate gives the
\code{ucpc.founder.proxy} quadratic \eqn{c'Kc}, corresponding to the
plug-in diversity proxy in Frontiers in Genetics equation (9).
It omits finite-panel contribution covariance and locus-specific shifts.
It is neither actual allele coancestry nor an input multiplicity for
\code{propsd}. Raw fractions outside \eqn{[0,1]} remain visible and make the
qualified approximation unavailable; they are never clipped. In-range
predictions can also be biased. An unresolved approximation cannot prove
infeasibility. Estimated numerical errors do not bound this biological
approximation bias. Without a cap these predictions do not change the search.

Optional scenarios generate actual immediate-DH marker-origin paths under
Haldane meiosis and select using
\eqn{I=\rho(G-\mu)/\sqrt{v}+\sqrt{1-\rho^2}\,\epsilon}, with independent
standard-normal noise and neutral ties. Direct genetic deviations preserve
ranking and response in the presence of a large common trait intercept.
The private common bank preserves candidate prefixes and leaves the caller's
R and simulation random streams unchanged. It runs at the final assessment,
not at each optimizer count candidate. No scenario score is inherited as a
prediction for changed counts.

With fixed reference frequencies \eqn{p_{0j}} and normalized weights,
\eqn{D=\sum_jw_jp_{0j}(1-p_{0j})}. Each selected replicate has
\deqn{q=\frac{\sum_jw_j(p_j-p_{0j})^2}{D},\qquad
H=2\sum_jw_jp_j(1-p_j).}
The \code{genomic.q} metric uses weighted VanRaden A/2 with all line pairs,
including the diagonal; \code{st.get.A} itself returns A. The selected
population is never used to recenter its reference. \code{gene.diversity}
is population allele diversity; individual DH heterozygosity is zero.
Nonlinear metrics are calculated within each replicate before averaging.
A zero reference denominator makes only the genomic quadratic unavailable.
Positive normalizers or marker weights that cannot be represented
consistently within the model's numerical precision produce a
\code{numeric.range} error.

Scenario variance summaries distinguish marginal selected-member variance,
expected finite-panel variance with denominator \eqn{R}, and expected sample
variance with denominator \eqn{R-1}. The marginal estimate is the mean
finite-panel variance plus the unbiased variance of replicate panel means;
it requires multiple completed replicates. Its SE is explicitly unavailable.
Other available SE values use independent whole-replicate summaries.
One retained individual has finite-panel variance zero and unavailable sample
variance; one replicate has unavailable SE. Empty selected populations have
reasoned unavailable selected summaries and no candidate draws. Resource
limits publish only whole completed replicates and identify a partial bank.
These fixed-model simulations describe the immediate DH generation; they
do not provide independent field validation.

}
\section{Additional DH cycle}{
The optional \code{genetic.model$future} group supports these named
components:
\describe{
\item{\code{policy}}{Mating policy \code{"random.disjoint"}.}
\item{\code{n.progeny}, \code{n.selected}}{Positive whole-number
produced and retained counts per next-cycle family, defaulting to 20 and
2, respectively.}
\item{\code{accuracy}}{Selection accuracy between zero and one,
defaulting to one.}
\item{\code{training}}{Named list of \code{replicates} and \code{seed},
defaulting to 32 and zero.}
\item{\code{evaluation}}{Named list of \code{replicates} and \code{seed},
defaulting to 256 and one.}
\item{\code{comparators}}{Optional named list of actual plans,
allocations, population evaluations or count tables; the default is
\code{NULL}.}
}
This group belongs to the additive inbred DH model described above. Current
production and retention counts come from the evaluated plan, allocation or
count table; joint planning budgets are described in
\code{\link{st.cross.progeny}}. The selected current
population contains an even number of at least two logical DH lines.
Uniform random disjoint matching uses every selected line exactly once.
Distinct lines with identical genotypes or from the same original family
may mate. New F1 crosses between those lines produce the next DH families.
A homozygous DH line has no segregation variance by itself.

The future family counts satisfy
\eqn{n_2\geq r_2\geq1}. Future selection accuracy is between zero and one.
The number of next-cycle families is half the current retained total;
their produced and retained totals are that number times \code{n.progeny}
and \code{n.selected}, respectively. Counts and replicate numbers are
supported positive whole numbers; seeds are supported nonnegative whole
numbers. Unknown, repeated and malformed fields produce errors.

The model uses the current frozen marker effects, actual parental
haplotypes, Haldane recombination and independent selection-index noise in
both cycles. Conditional next-cross variance is calculated from the
actual selected genotypes. Mean conditional variance and mean conditional
SD are separate summaries; neither a kinship matrix nor the square root
of mean variance replaces those genotypes.

\code{criterion="future"} requires joint planning, an explicit
\code{progeny} group and a valid requested genetic model. Its objective is
the mean genetic value in the next selected population. Search compares
only complete evaluations on a fixed training bank. A separate evaluation
bank assesses the chosen plan. Training and evaluation streams differ even
when their user seeds are equal. A finite-bank optimum is not an optimum
of the unknown expectation. Larger searches retain their best-found status
and unresolved gap. Existing current-generation feasibility restrictions,
including an active selected-founder proxy cap, still apply to actual
current counts. Such a proxy cap is not a bound on actual next-generation
genomic coancestry.

An explicit \code{future} group with an ordinary mean or usefulness plan,
\code{st.cross.evaluate}, or \code{st.cross.allocate} requests a final
assessment only. It does not alter that operation's current-generation
objective or search path. Odd or insufficient current retention yields an
unavailable future assessment without changing current feasibility;
future optimization rejects those retention totals before search.

\code{comparators} is an optional named list of actual plans, allocations,
population evaluations or count tables. Names are unique and nonempty; \code{selected} and
\code{G1.reference} are reserved. Comparators are checked against the same
genetic core, current resources and applicable current restrictions.
Unavailable or invalid comparator status is individual. It cannot change
the chosen plan's feasibility. Future planning also retains its checked
current-generation reference. Chosen-minus-comparator differences and SE
use whole shared evaluation replicates. Evaluation outcomes do not rank
search proposals or select later warm starts.

Normalized model metadata retains only the future mating policy, family
counts and accuracy. A new explicit request is required to execute an
assessment. Banks, comparator commands, results and certificates from a source
object are not inherited. The frozen genetic snapshot and its signature are unchanged
by future policy metadata. Planner starts contribute counts and pairs,
without supplying an implicit model or future policy.

One cumulative genetic cell/work budget covers preparation, caches,
training and assessment, including all frontier points. A finite
assessment allowance is reserved before search; unusual costs can still
leave it incomplete. Requested and completed whole replicates and the
resource status remain explicit. Partial training banks do not compete
with complete candidates. Genomic coancestry, gene diversity and variance
summaries use the same fixed-reference definitions as the current DH
assessment, with nonlinear statistics calculated per replicate.

The forecast is conditional on fitted effects, the declared selection
model and one further cycle. Simulation agreement is model consistency,
not independent field validation or an indefinite breeding guarantee.
}
\value{
An object of class \code{st.cross.population}, with normalized \code{crosses},
\code{population} summaries, \code{contributions}, \code{relatedness},
\code{feasibility}, \code{availability}, \code{diagnostics} and declared
\code{inputs}. Original input-row mapping and source status remain available.
The print method returns its input invisibly.

The \code{future} group is a named input rather than a separate R function.
Activated planning and assessment results contain a \code{future} group
with the horizon, policy and resources, training and evaluation identities,
completion status, current and next selected-population summaries,
conditional variance/SD, paired comparator differences and resource
diagnostics. Their exact availability depends on completed evaluations.
Existing current-generation analytic fields retain their meanings.
For an explicit assessment, \code{population$future.mean} contains the
same available held-out G2 mean, and \code{availability$G2} and its reason
describe its availability and completion. Uncertainty remains in the
dedicated \code{future} summaries.
}
\references{
Okoye UJ, Frisch M, Herzog E (2026). Optimization of family sizes in sets
of crosses with a greedy allocation strategy based on automatic
differentiation. Frontiers in Plant Science 17:1727383.
\doi{10.3389/fpls.2026.1727383}.

Allier A, Moreau L, Charcosset A, Teyssedre S, Lehermeier C (2019).
Usefulness Criterion and Post-selection Parental Contributions in
Multi-parental Crosses: Application to Polygenic Trait Introgression.
G3: Genes, Genomes, Genetics 9:1469--1479.
\doi{10.1534/g3.119.400129}.

Allier A, Lehermeier C, Charcosset A, Moreau L, Teyssedre S (2019).
Improving Short- and Long-Term Genetic Gain by Accounting for Within-Family
Variance in Optimal Cross-Selection. Frontiers in Genetics 10:1006.
\doi{10.3389/fgene.2019.01006}.

Niehoff TAM, ten Napel J, Bijma P, Pook T, Wientjes YCJ, Hegedűs B,
Calus MPL (2024). Improving selection decisions with mating information by
accounting for Mendelian sampling variances looking two generations ahead.
Genetics Selection Evolution 56:41.
\doi{10.1186/s12711-024-00899-2}. Equations (2), (8), and (9) motivate
the conditional-variance distinction; the DH cycle here uses actual
selected genotypes.

Hamazaki K, Iwata H (2024). AI-assisted selection of mating pairs through
simulation-based optimized progeny allocation strategies in plant breeding.
Frontiers in Plant Science 15:1361894.
\doi{10.3389/fpls.2024.1361894}. Sections 1 and 2.1 concern optimization
of future progeny populations.

Niehoff TAM, ten Napel J, Calus MPL (2025). Practical Considerations When
Using Mendelian Sampling Variances for Selection Decisions in Genomic
Selection Programs. Journal of Animal Breeding and Genetics 142:419--437.
\doi{10.1111/jbg.12913}. Sections 4.1 and 4.4 discuss accuracy and
model-explainable variance under estimated marker effects.
}
\seealso{\code{\link{st.cross.allocate}}, \code{\link{st.cross.plan}},
\code{\link{st.cross.progeny}}, \code{\link{propsd}},
\code{\link{st.cross.frontier}}, \code{\link{st.cross.future}}
}
