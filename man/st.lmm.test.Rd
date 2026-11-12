\name{st.lsm-hypotheses}
\alias{st.lsm-hypotheses}
\title{Fixed-Effect Hypotheses Through st.lsm}
\description{The direct coefficient path in \code{st.lsm} tests
\eqn{L\beta=\mathrm{rhs}} for supported moment, ML and REML fits.}
\details{
\code{L} accepts a vector or matrix; \code{joint=TRUE} requests a joint
F test. Named columns align with the fixed coefficients. A scalar or
per-row \code{rhs} specifies nonzero null values. Consistent redundant
rows are reduced; non-estimable or inconsistent hypotheses are rejected.
Kenward--Roger uses adjusted REML covariance, hypothesis-specific denominator
DF and F scaling. Satterthwaite combines orthogonalized one-dimensional
moment approximations. \code{ddf="asymptotic"} selects a Wald F test with
infinite denominator DF; \code{test="Chisq"} selects the global asymptotic
Wald chi-squared test.
Without \code{joint=TRUE}, the result stores the linear estimates and
covariance and supports ordinary summaries and simultaneous inference.
Symbolic linear equalities and inequalities are also accepted in \code{L}.
The existing least-squares mean construction is unchanged; an existing means or
contrast object can be passed to \code{st.lsm(e, joint=TRUE)} for the full
model-based joint calculation.
}
\value{Both direct joint tests and row-wise construction return an
\code{st.lsm} object. For a direct joint F test, \code{summary} returns the
stored table with \code{DF1}, \code{DF2}, \code{F}, \code{F.scaling},
\code{p} and \code{method}.}
\seealso{\code{\link{st.lsm}}, \code{\link{st.lmm}}}
