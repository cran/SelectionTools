\name{st.ocs}
\alias{st.ocs}
\alias{print.st.ocs}
\title{Optimum Parental Contributions}
\description{
Finds parental gene contributions that maximize supplied additive merit under
an upper limit on group relatedness and lower and upper contribution bounds.
}
\usage{
st.ocs(merit, K, theta, lower = 0, upper = 1, control = list())

\method{print}{st.ocs}(x, digits = NULL, top = NULL, bottom = NULL,
                     char.per.line = NULL, ...)
}
\arguments{
  \item{merit}{Finite numeric vector of parental additive merits. Each entry
  corresponds to the same position in both dimensions of \code{K}.}
  \item{K}{Finite symmetric positive-semidefinite relationship matrix on the
  scale used by \code{theta}. For an additive relationship matrix \code{A},
  the kinship convention is \code{K = A/2}. Singular matrices are supported.
  Matrix names are either absent or supplied for both dimensions and equal.
  If \code{merit} is also named, its names must match these names in the
  same order.}
  \item{theta}{Finite upper limit on \eqn{c'Kc}, where \eqn{c} is the vector
  of parental gene contributions. The meaning and reference population of
  this limit follow from the supplied matrix.}
  \item{lower, upper}{Lower and upper contribution bounds, each supplied as
  an unnamed numeric scalar or as a vector of the same length as \code{merit}.
  Bounds satisfy \code{0 <= lower <= upper <= 1}. Full-length named bounds
  follow the parental order.}
  \item{control}{Named list of numerical controls. \code{tolerance}
  (default \code{1e-8}) controls numerical convergence; \code{maxit}
  (default \code{10000}) limits inner optimization iterations;
  \code{max.outer} (default \code{100}) limits iterations per outer search
  phase; \code{rank.tolerance} (default \code{1e-10}) controls numerical
  rank decisions. The method scales the merit vector and matrix internally;
  convergence tolerances apply on those scales. Diagnostic residuals and
  objective bounds are returned in the original input units. Unknown controls
  produce an error.}
  \item{x}{An optimum-contribution result.}
  \item{digits, top, bottom, char.per.line}{Optional printing controls.
  Unspecified values follow the corresponding \code{st.options} settings.}
  \item{...}{Additional printing arguments.}
}
\details{
The optimization problem is
\deqn{
\mathop{\mathrm{maximize}}_c\quad g'c
\quad\mathrm{subject\ to}\quad
\mathbf{1}'c=1,\quad l\leq c\leq u,\quad c'Kc\leq\theta,
}{
maximize g'c subject to sum(c)=1, lower<=c<=upper, c'Kc<=theta,
}
where \eqn{g} is \code{merit}. The returned \code{gain} is the weighted
mean of the supplied merits.

Feasibility is checked using the contribution bounds and the minimum
attainable relatedness. An infeasible request is reported explicitly;
constraints are not relaxed. The numerical method also handles tied merits,
active bounds and singular positive-semidefinite matrices. An iteration limit
has a separate status from a converged solution.

The matrix scale remains fixed across solutions. A centered genomic matrix
defines a relatedness index relative to its allele-frequency reference.
Dividing an additive relationship matrix by two does not by itself turn this
index into an absolute probability of identity by descent. Negative
off-diagonal entries can occur in a centered genomic matrix.

The contributions are nominal parental gene shares before within-family
selection. They do not specify parent pairs, integer family sizes or realized
ancestry after selection. In a non-selfed crossing plan, a contribution above
one half cannot be realized; even contributions bounded by one half may need
additional pair and count feasibility checks. Such mating restrictions are
separate from the contribution problem solved here.

Input vectors and matrices and the caller's thread settings are preserved.
Calculations use the supplied numerical precision without rounding results
through formatted text. Numerical tolerances and certificates describe
floating-point convergence.
}
\value{
A list of class \code{st.ocs} with the following components. Its printed
representation summarizes the fit and the largest and smallest parental
contributions under the current printing options.
\item{contributions}{Named numeric vector in parental order, or \code{NULL}
when the request is infeasible or no feasible incumbent has been found within
the numerical limits.}
\item{gain}{Weighted parental merit, \eqn{g'c}, for the returned solution.}
\item{relatedness}{Achieved group relatedness, \eqn{c'Kc}.}
\item{minimum.relatedness}{Relatedness from the minimum-relatedness search
under the contribution bounds. The difference from
\code{diagnostics$minq.lower.bound} measures the remaining numerical gap in
that search. The main status describes the merit optimization; a certified
merit optimum can be available before this separate gap has converged.}
\item{theta}{Requested relatedness limit.}
\item{status}{\code{"optimal"} for a certified ordinary solution,
\code{"minimum.relatedness"} for a certified boundary solution,
\code{"infeasible"} for an infeasible request, or \code{"maxit"} when
the iteration or numerical-certification limits prevent convergence.}
\item{diagnostics}{Numerical certificates, constraint residuals, objective
bounds, iteration counts and matrix-rank information. A minimum-relatedness
boundary may have a separate certificate because a finite relatedness
multiplier need not exist there.}
\item{bounds}{Expanded lower and upper contribution bounds.}
\item{control}{Numerical control values used for the fit.}
}
\references{
Wellmann R (2019) Optimum contribution selection for animal breeding and
conservation: the R package optiSel. BMC Bioinformatics 20:25.
}
\seealso{
\code{gs.esteff.rr}, \code{gs.estimate.gv}, \code{st.get.A},
\code{gs.cross.info}
}