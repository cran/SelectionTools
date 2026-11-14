\name{st.lmm.diagnostics}
\alias{st.lmm.diagnostics}
\alias{hatvalues.st.lmm}
\alias{rstandard.st.lmm}
\title{Leverage and Residual Diagnostics for Gaussian Mixed Models}
\description{
Calculates marginal or conditional residuals, residual variances,
standardized residuals and leverage at the fitted covariance parameters.
The diagnostic table retains explicit response-cell, row and trait keys.
}
\usage{
st.lmm.diagnostics(object, type = "conditional")
\method{hatvalues}{st.lmm}(model, type = "conditional", ...)
\method{rstandard}{st.lmm}(model, type = "conditional", ...)
}
\arguments{
\item{object,model}{An \code{st.lmm} object fitted by ML or REML.
Univariate and multivariate responses are supported. The \code{hatvalues}
and \code{rstandard} methods also support univariate method-of-moments
models containing only fixed effects.}
\item{type}{\code{"conditional"} includes fitted random effects in the
mean. \code{"marginal"} uses the fixed mean, including fitted offsets.}
\item{...}{Additional arguments to the generic; currently unused.}
}
\details{
All quantities condition on the fitted covariance parameters. They do not
add covariance-parameter estimation uncertainty or apply Kenward--Roger
adjustments. The result records whether the original covariance fit
converged; evaluating diagnostics does not make an unconverged fit valid.

The covariance of the observed response cells is \eqn{V=G+R}, where
\eqn{G} is the contribution of all random effects and \eqn{R} is the
residual covariance. The matrix \eqn{X} contains an independent basis for
the fixed mean, with \eqn{C=(X'V^{-1}X)^{-1}} and
\eqn{P=V^{-1}-V^{-1}XCX'V^{-1}}.

Marginal leverage is the diagonal of \eqn{H_m=XCX'V^{-1}}.
Conditional leverage is the diagonal of \eqn{H_c=I-RP}.
With correlated residuals these matrices need not be symmetric; leverage
values are not truncated to the interval from zero to one.

The conditional residual is the observed response minus its fitted fixed
and random mean, with variance given by the diagonal of \eqn{RPR}.
The marginal residual subtracts only the fixed mean and offset, with
variance given by the diagonal of \eqn{V-XCX'}.
\code{rstandard} divides the selected residual by the square root of its
variance. A zero variance gives \code{NA} and an explicit status rather
than a zero standardized residual.

For a univariate fixed-only method-of-moments model, \code{hatvalues}
uses the weighted fixed-design leverage and \code{rstandard} uses the
ordinary regression residual variance with DF equal to the number of
observations minus the fixed-design rank. Positive precision weights
multiply the response residual by their square root before
standardization. Marginal and conditional results agree in this case.
The \code{st.lmm.diagnostics} table remains available only for ML and REML.

\code{residuals(object)} uses the conditional response residual for a
likelihood mixed model and the ANOVA projection residual for a mixed
moments fit. The corresponding explicit types are \code{"conditional"}
and \code{"projection"}. For a model without random effects, the default
is the residual from the fixed mean.
\code{residuals(object, type="marginal")} subtracts the fixed mean.
For ML and REML, \code{type="pearson"} divides the conditional response
residual by the square root of the diagonal of \eqn{R}. Pearson residuals
and \code{rstandard} generally use different denominators.

\code{plot(object, which=1:4)} displays residuals versus fitted values,
a normal Q--Q plot of standardized residuals, a scale-location plot, and
standardized residuals versus leverage. Character selectors are
\code{"residuals"}, \code{"qq"}, \code{"scale-location"} and
\code{"leverage"}. A model with subset-selection results retains its
subset inclusion plot by default; \code{which} requests a
diagnostic panel. Plots are descriptive, since residuals can be correlated.
They report unavailable cells, annotate at most \code{id.n} observations,
and invisibly return their exact numeric plot data. Plotting does not change
the graphics layout; the usual graphics controls determine a multi-panel page.
}
\value{
\code{st.lmm.diagnostics} returns a data frame with columns
\code{cell}, \code{row}, \code{original.row}, \code{row.name},
\code{trait}, \code{observed}, \code{estimable}, \code{response},
\code{fitted}, \code{residual}, \code{variance}, \code{standardized},
\code{leverage} and \code{status}. It describes retained response cells,
including missing responses. The original-row key refers to the row of the
data supplied to the fitted model; it can repeat when the fitted subset
contains repeated rows. Attributes state the diagnostic type, covariance
assumption, leverage and variance definitions, and convergence status.

\code{hatvalues} and \code{rstandard} return matrices with one row per
observation and one column per trait. They preserve row and
trait labels, and restore omitted selected rows with \code{NA} when the
model uses \code{na.exclude}. Undefined or nonestimable quantities remain
\code{NA}. The diagnostic table itself uses retained-row keys and is not
expanded by \code{na.exclude}.
}
\seealso{
\code{\link{st.lmm}}, \code{\link{st.lmm.components}},
\code{\link{influence.st.lmm}}, \code{\link{rstudent.st.lmm}}
}
\keyword{models}

