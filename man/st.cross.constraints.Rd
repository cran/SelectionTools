\name{st.cross.constraints}
\alias{st.cross.constraints}
\title{Practical Restrictions for Cross Planning and Progeny Allocation}
\description{
The \code{constraints} argument of \code{st.cross.plan},
\code{st.cross.allocate} and \code{st.cross.evaluate} describes parent
counts, reproductive roles, capacities and nominal population restrictions.
Every request names its population basis and intersects existing restrictions.
}
\details{
The optional named list has five supported components:
\describe{
  \item{\code{parents}}{A table of unique character \code{parent} IDs,
  optional character \code{pool}, and optional logical \code{seed} and
  \code{pollen} eligibility. An omitted eligibility column means that role
  is allowed. A supplied table covers the candidate universe; reference-only
  matrix names do not create additional candidates. One pool per parent is
  supported. Pool labels alone do not declare reproductive direction.}
  \item{\code{limits}}{A table with \code{parent}, \code{pool},
  \code{role}, \code{basis}, \code{measure}, \code{lower} and \code{upper}.
  A missing parent selects a group; a missing pool selects all pools.
  Roles are \code{"all"}, \code{"seed"} and \code{"pollen"}.
  Measures are \code{"parents"} for distinct used parents,
  \code{"uses"} for integer incidences, and \code{"share"} for normalized
  contribution. Parent counts apply to groups. Equal lower/upper bounds
  specify exact limits. Integer bounds use supported whole numbers;
  shares lie in \eqn{[0,1]}. Repeated identical selectors are errors.}
  \item{\code{equal.use}}{A table of \code{pool}, \code{role} and
  \code{basis}. All actually used parents in each selected subset have
  equal integer incidence. It does not require every eligible parent to
  be used; a distinct-parent limit specifies that requirement separately.}
  \item{\code{relatedness}}{A table of \code{basis} and \code{max},
  imposing \eqn{c'Kc <= \theta} on overall nominal shares with the original
  named \code{K}. It does not replace the existing \code{theta} argument.}
  \item{\code{diversity}}{A list with \code{G}, \code{max} and
  \code{sample}. Samples are \code{"family.slots"},
  \code{"produced.founder.slots"} and
  \code{"retained.founder.slots"}. Literal integer founder multiplicities
  and the complete fixed covariance reference are used by \code{propsd}.
  Its estimated proportional loss is strictly below \code{max}.
  This request intersects the existing family-slot diversity restriction.}
}

The supported incidence bases are \code{"family.slots"},
\code{"produced"}, \code{"retained"} and \code{"operations"}.
For produced, retained or identified operation weights \eqn{w_a}, with
\eqn{T=\sum_a w_a}, seed and pollen incidences are
\deqn{S_i=\sum_{a:\mathrm{seed}(a)=i}w_a,\qquad
L_i=\sum_{a:\mathrm{pollen}(a)=i}w_a.}
Overall gene share is \eqn{(S_i+L_i)/(2T)}; seed and pollen shares are
\eqn{S_i/T} and \eqn{L_i/T}. A pool sums these quantities with the same
overall or role denominator; it is not normalized within its own pool.
An overall distinct-parent count counts a parent once even when it occurs
in both roles. Role-specific counts count its use in each selected role.

Family slots count each active unordered biological family once, regardless
of duplicated or reciprocal rows. Its overall denominator is twice the
number of active families. Role-specific family slots require an unambiguous
direction for every active family. Multiple reproductive records can contribute
to one pooled biological response, but cannot create additional families or
selection response by duplication.

Generated role-aware plans and allocations carry a \code{directions} table
with \code{family}, \code{seed} and \code{pollen}, without duplicated
progeny counts. Its actual unordered parent key joins to the current
\code{crosses} counts. An old role-free object's display order does not
identify reproductive direction. A raw data frame with an explicit new
role contract declares \code{P1Name} as seed and \code{P2Name} as pollen.
Its original reciprocal rows retain their separate actual counts for
reproductive accounting; response pools by biological family.

Operation multiplicities are identified only by explicit
\code{n.operations}; they are not inferred from progeny counts. Each is a
nonnegative exact whole number at most \eqn{2^{53}-1}, and their total obeys
the same limit. A positive produced count requires positive identified
operation count when that basis is requested. Zero production with recorded
positive operations describes an unsuccessful operation. The current planner
and fixed-family allocator do not optimize operation multiplicities or
reciprocal production splits.

Produced and retained counts keep the selection kernel's per-family integer
domain. Group incidence can exceed \eqn{2^{53}-1}; feasibility and exact
character reports retain its integer value. Numeric report fields are
unavailable when the value cannot be represented in the supported exact
numeric domain. Share checks use exact incidence ratios and the represented
numeric bounds. Integer equal-use restrictions describe exact equal use.
The literal multiplicity total for a PropSD sample is separately limited
to \eqn{2^{53}-1}; a larger sample is not silently rescaled.

No-selection and zero-total situations retain their declared meanings.
Shares or relatedness with a zero denominator are unavailable. The
equal-family planner has no retained or operation counts to infer, so
restrictions requiring those counts are rejected. Fixed-family allocation
with produced-share, capacity, relatedness or produced-diversity coupling
uses a constrained search. The separable greedy/dynamic-programming
optimality claim does not transfer to that problem. Search limits and
unresolved numerical boundaries remain explicit.

Nominal founder-slot PropSD and \eqn{c'Kc} are not actual selected-DH
genomic diversity.
}
\value{
This topic describes an argument rather than a separate function.
Opted-in results retain effective practical restrictions under
\code{constraints$practical}; the current explicit request remains under
\code{inputs$practical}. Recomputed checks are in
\code{diagnostics$practical}. Evaluation reports independently recomputed
feasibility. Reports distinguish passed, failed, numerically
unresolved and unavailable restrictions; a saved source status does not
establish feasibility or optimality.
}
\references{
Gorjanc G, Hickey JM (2018). AlphaMate: a program for optimizing selection,
maintenance of diversity and mate allocation in breeding programs.
Bioinformatics 34:3408--3411. \doi{10.1093/bioinformatics/bty375}.
The associated manual's named mating/crossing-constraint sections describe
parent counts, role-specific counts, equal contributions and contribution
bounds: \url{https://github.com/AlphaGenes/AlphaMate/blob/main/doc/AlphaMate.txt}.
}
\seealso{\code{\link{st.cross.plan}}, \code{\link{st.cross.allocate}},
\code{\link{st.cross.evaluate}}, \code{\link{propsd}}}
