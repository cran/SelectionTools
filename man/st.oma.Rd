\name{st.oma}
\alias{st.oma}
\alias{print.st.oma}
\title{OMA: Optimum Mate Allocation}
\description{
Chooses continuous allocations to allowed diploid matings to maximize their
weighted merit under parental and mating bounds and two-generation relatedness
constraints. OMA denotes optimum mate allocation.
}
\usage{
st.oma(parents, matings, K, dF, base = "current", control = list())

\method{print}{st.oma}(x, digits = NULL, top = NULL, bottom = NULL,
                     char.per.line = NULL, ...)
}
\arguments{
  \item{parents}{Data frame containing parental IDs and lower and upper
  contribution bounds. Its three columns are \code{id}, \code{min},
  and \code{max}, in that order. IDs are unique, complete character strings. Bounds are
  finite and satisfy \code{0 <= min <= max <= 1}.}
  \item{matings}{Data frame containing allowed pairs, their fixed merit,
  and lower and upper allocation bounds. Its five columns are
  \code{parent1}, \code{parent2}, \code{merit}, \code{min}, and
  \code{max}, in that order. Parent IDs match \code{parents} and \code{K}.
  Self matings are supported. Each unordered pair occurs once; reversed
  duplicate rows produce an error.
  Merits are finite numeric values; bounds satisfy
  \code{0 <= min <= max <= 1}.}
  \item{K}{Finite symmetric positive-semidefinite diploid kinship matrix.
  Both dimensions have the same unique, complete parental names in the
  same order. It may contain additional reference individuals that do not
  occur in the allowed matings.}
  \item{dF}{Finite numeric vector of length one or two specifying a scalar
  change or lower and upper changes relative to the reference value.
  Both endpoints are strictly less than one; negative endpoints are supported.
  Two values are ordered from lower to upper. One value supplies both
  endpoints, giving an equality for the first-generation constraint.}
  \item{base}{Either \code{"current"} or \code{"RM"}. The former
  uses average current inbreeding as the reference; the latter uses a
  randomly mated reference including selfing. The supplied matrix and
  reference remain fixed throughout optimization.}
  \item{control}{Named list of numerical controls. Unknown and repeated
  names produce errors. The supported components are:
  \describe{
    \item{\code{maxit}}{Positive whole-number iteration limit no larger
    than \code{.Machine$integer.max};
    default \code{500}.}
    \item{\code{tolerance}}{Finite numerical tolerance for the reported
    scaled primal residual and normalized merit-bound gap, strictly between
    zero and one;
    default \code{1e-8}.}
  }}
  \item{x}{An optimum mate-allocation result.}
  \item{digits, top, bottom, char.per.line}{Optional printing controls.
  Unspecified values follow the corresponding \code{st.options} settings.}
  \item{...}{Additional printing arguments.}
}
\details{
For mating \eqn{e=(i,j)}, the continuous allocation \eqn{x_e} is the fraction
of produced progeny assigned to that mating. Allocations sum to one. A
parent's gene contribution is
\deqn{
c_i=\frac{1}{2}\sum_e x_e
\left[\mathbf{1}(i_e=i)+\mathbf{1}(j_e=i)\right].
}{
c[i] = sum of allocations across both parental slots occupied by i divided by two.
}
A self mating contributes both parental slots to the same parent.
The objective is the allocation-weighted merit \eqn{h'x}. A supplied merit
can represent a predicted cross mean or a fixed usefulness score. Score
calculation and the biological selection scenario are determined before
allocation; changing allocation alone does not recalculate the score.
A subsequent resource allocation with different family sizes requires scores
for those sizes. Contribution constraints concern nominal gene shares before
within-family selection, rather than realized selected-descendant ancestry.

For diploids, the first- and second-generation quantities are
\deqn{F_1=\sum_e K_{i_ej_e}x_e,\qquad F_2=c'Kc.}{
F1 = sum(K[parent1,parent2] * allocation); F2 = sum over i,j of c[i]*K[i,j]*c[j].
}
With \code{base = "current"}, the reference is
\eqn{F_0=\operatorname{mean}(2\operatorname{diag}(K)-1)}.
With \code{base = "RM"}, it is \eqn{F_0=\operatorname{mean}(K)}.
These means use the complete supplied matrix before restriction to eligible
mating parents. Extra reference individuals can therefore change the
constraints even when they do not occur in any selected mating.

For \code{dF = c(dL, dU)}, the constraints are
\deqn{
d_L+(1-d_L)F_0\leq F_1\leq d_U+(1-d_U)F_0,
}{
dL + (1 - dL) * F0 <= F1 <= dU + (1 - dU) * F0,
}
\deqn{
F_2\leq d_U(2-d_U)+(1-d_U)^2F_0.
}{
F2 <= dU * (2 - dU) + (1 - dU)^2 * F0.
}
The reported transformed quantities are
\eqn{dF_1=(F_1-F_0)/(1-F_0)} and
\eqn{dF_2=1-\sqrt{(1-F_2)/(1-F_0)}} where these are defined.
A scalar \code{dF} is an equality for \eqn{F_1}, rather than an upper
bound alone. The current reference is undefined when \eqn{F_0=1}; the
randomly mated reference is the relevant alternative for fully inbred lines.

A centered genomic matrix describes a relatedness index relative to its
allele-frequency reference. The formulas do not by themselves make that
index an absolute identity-by-descent probability or an observed biological
inbreeding rate. The first-generation average within-mating quantity and
second-generation group quantity are distinct constraints.

Continuous mate allocation does not supply a fixed number of distinct
families or integer progeny counts. Optimum cross selection with equal
positive integer family sizes is provided by \code{st.cross.plan}.
Conversion of continuous allocations into an operational integer plan is
a separate constrained problem; allocations are not rounded automatically.
}
\value{
An object of class \code{st.oma} with the following components:
\describe{
  \item{allocations}{Numeric vector of continuous allocation fractions in
  the input mating order. The values sum to one for a feasible result.
  No output thresholding or integer rounding is applied.}
  \item{contributions}{Named vector of parental gene contributions in
  the input parental order.}
  \item{parents}{The supplied parental table with its contribution column.}
  \item{matings}{The supplied mating table with allocation and kinship
  columns. All candidate rows are retained, including zero allocations.}
  \item{F0, F1, F2}{The full-reference value, allocation-weighted mating
  kinship and parental group quantity, respectively.}
  \item{dF1, dF2}{The corresponding transformed first- and second-generation
  changes, where defined.}
  \item{merit}{The allocation-weighted supplied cross merit.}
  \item{status}{\code{"optimal"} when the reported scaled primal residual
  and normalized merit-bound gap satisfy the requested tolerance. The gap
  is divided by the largest absolute deviation of supplied merits from
  the midpoint of their minimum and maximum, with a divisor of one for
  constant merits. The reported merit and gap retain their original units.
  Other statuses are \code{"infeasible"} when infeasibility is established, or
  \code{"maxit"} when the iteration or numerical limit leaves the
  request unresolved. A \code{"maxit"} result can retain a feasible
  incumbent. Without a feasible incumbent, allocations, contributions and
  their derived values are missing rather than a feasible mating plan.}
  \item{diagnostics}{Numerical reason, iteration count, constraint
  residuals, infeasibility lower bound, matrix eigenvalue and curvature
  information, reference size and problem dimensions. The overall primal
  residual combines absolute allocation/simplex/contribution residuals
  with scaled first- and second-generation violations. Individual
  \code{F1.residual} and \code{F2.residual} values use the original
  kinship units. The \code{affine.roundoff} component records numerical
  uncertainty in the scaled first-generation equality. Components
  \code{last.candidate.residual}, \code{last.candidate.F1.residual},
  and \code{last.candidate.F2.residual} describe the latest candidate
  assessed, which need not be a feasible returned incumbent.}
  \item{bounds}{The requested two-element \code{dF} interval, transformed
  \code{F1.lower}, \code{F1.upper} and \code{F2.upper}, merit upper bound,
  remaining gap and bound description. The gap is evaluated using centered
  merits before the common merit offset is restored. This preserves the
  same simplex objective while avoiding cancellation under a large common
  offset. The original-unit \code{merit.upper} remains rounded upward;
  subtracting the rounded reported merit from it can therefore give a
  larger difference than the reported gap.}
  \item{control}{The iteration limit and tolerance used.}
  \item{inputs}{The supplied parent and mating tables, complete matrix,
  \code{dF} request and reference convention retained for reconstruction.}
}
The print method returns \code{x} invisibly.
}
\references{
Endelman JB (2025) Genomic prediction of heterosis, inbreeding control,
and mate allocation in outbred diploid and tetraploid populations.
Genetics 229:iyae193. \doi{10.1093/genetics/iyae193}.
}
\seealso{\code{\link{st.ocs}}, \code{\link{st.cross.plan}},
\code{\link{gs.cross.info}}, \code{\link{st.options}}}
