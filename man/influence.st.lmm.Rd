\name{influence.st.lmm}
\alias{influence.st.lmm}
\alias{rstudent.st.lmm}
\alias{cooks.distance.st.lmm}
\alias{dfbeta.st.lmm}
\alias{dfbetas.st.lmm}
\alias{summary.st.lmm.influence}
\alias{print.st.lmm.influence}
\alias{plot.st.lmm.influence}
\alias{cooks.distance.st.lmm.influence}
\alias{dfbeta.st.lmm.influence}
\alias{dfbetas.st.lmm.influence}
\title{Exact Deletion Influence and Predictive Studentization}
\description{
Refits a Gaussian mixed model after deleting cases or groups, preserving
the calibrated fixed and random designs and covariance levels. Reports
coefficient and covariance changes, Cook's distance, DFBETAs, predictive
residuals, and the status of every requested deletion.
}
\usage{
\method{influence}{st.lmm}(model, group = NULL, units = NULL,
  control = list(), keep.fits = FALSE, ...)
\method{rstudent}{st.lmm}(model, influence = NULL,
  control = list(), ...)
\method{cooks.distance}{st.lmm}(model, influence = NULL, ...)
\method{dfbeta}{st.lmm}(model, influence = NULL, ...)
\method{dfbetas}{st.lmm}(model, influence = NULL, ...)
\method{summary}{st.lmm.influence}(object, ...)
\method{print}{st.lmm.influence}(x, ...)
\method{plot}{st.lmm.influence}(x, which = c("cook", "dfbetas"),
  parm = NULL, labels = NULL, id.n = 3L, main = NULL, col = NULL, ...)
\method{cooks.distance}{st.lmm.influence}(model, ...)
\method{dfbeta}{st.lmm.influence}(model, ...)
\method{dfbetas}{st.lmm.influence}(model, ...)
}
\arguments{
\item{model}{For fitted-model methods, a converged \code{st.lmm} object
fitted by ML or REML. For influence-object methods, an
\code{st.lmm.influence} object.}
\item{group}{NULL deletes one retained data row per case, including all
its observed traits. A single character string names a column of the
retained data whose values define groups. Otherwise supply a nonmissing
atomic vector with exactly one value per retained row. A vector indexed
to the original input rows is not accepted when its length differs from
the retained data.}
\item{units}{NULL requests every deletion unit. Otherwise supply distinct
unit IDs or labels from the inventory. IDs follow the first appearance
of each unit in the retained data. Units not requested remain visible
in the result.}
\item{control}{Named fitting controls that override the original controls
for deletion fits. The original ML or REML method is retained.}
\item{keep.fits}{Retain complete deletion models when TRUE. Raw coefficient
and covariance estimates, diagnostic status and warnings are retained
regardless of this choice when a fit was returned.}

\item{influence}{NULL computes deletion influence. Alternatively supply
an influence result from the same fitted model. \code{rstudent} requires
a complete case-deletion result; grouped or partial results are rejected.}
\item{object,x}{An \code{st.lmm.influence} object.}
\item{which}{Influence panels: \code{"cook"} and \code{"dfbetas"}, or
their integer indices 1 and 2. Several distinct panels may be requested.}
\item{parm}{NULL plots the maximum absolute DFBETA per deletion unit.
Fixed-coefficient names or indices request a separate signed DFBETA
panel for each selected coefficient.}
\item{labels}{Optional character labels, one per inventory unit.
The default uses the unit labels. Missing labels are not annotated.}
\item{id.n}{Maximum number of units to annotate in each panel, selected
by the largest absolute displayed influence. Zero suppresses annotations.}
\item{main}{One title, or one title per requested panel.}
\item{col}{One colour, or one colour per inventory unit.}
\item{...}{For fitted-model \code{cooks.distance}, \code{dfbeta} and
\code{dfbetas} without a supplied influence object, arguments to
\code{influence}. For plotting, graphical parameters such as
\code{xlim}, \code{ylim}, \code{xlab}, \code{ylab}, \code{pch} and
\code{cex}. For printing, arguments to the result-table printer.
Other methods currently ignore additional arguments.}
}
\details{
Repeated fits run in the current R process. Numerical kernels use the
existing OpenMP and BLAS thread controls; \code{\link{st.set.num.threads}} describes these controls.
Intermediate state stays in C-owned memory.

\subsection{Deletion fitting}{
Deletion marks the requested observed response cells missing in the
compiled model. The full fixed design, random and residual coefficient
bases, factor levels, known relationships and covariance kernels remain
calibrated to the original fit. Predictor transformations are not reevaluated.
Every free covariance parameter is re-estimated by the original ML or
REML method. Fixed covariance blocks retain their supplied values.

A disappearing grouping level remains available for prediction, including
its relationships to observed levels in a known covariance matrix.
Deleting all observed support for a free covariance stratum instead gives
an explicit \code{"lost covariance stratum"} result. The function does
not silently remove or fix that covariance parameter. With \code{homdiag}
or homogeneous compound-symmetric or AR1 residual covariance, an occasion
may lose all its observations if other occasions still inform the common
scale. The fitted occasion positions and AR1 gaps are retained. Losing all
support for a separate occasion variance or an entire common-scale stratum
remains a lost-stratum failure.

Every requested unit retains its warnings, convergence and rank information.
Failed fits, nonconvergence, nonestimable held-out means, changed fixed rank,
zero prediction-error variance and insufficient observations are reported
explicitly. A failure does not stop unrelated deletion units. Failed or
undefined inferential diagnostics remain \code{NA}; raw estimates can still
be inspected when the fitting engine returned them.
}

\subsection{Coefficient changes}{
Changes use the original estimate minus the deletion estimate.
\code{dfbeta} returns fixed-coefficient changes, only for individual
coefficients estimable in both fits. \code{dfbetas} divides each such
change by that coefficient's deleted-fit standard error. These standard
errors use the unadjusted fitted covariance, even if the original model
requested Kenward--Roger inference.

Cook's distance is the squared change in the estimable fixed-coefficient
vector, measured with the inverse original unadjusted fixed-coefficient
covariance, divided by the original fixed rank. Deleted coefficients are
mapped into the original independent fixed basis before this calculation.
Cook's distance is undefined if the estimable fixed-effect space changes,
the deletion does not converge, or the original fixed rank is zero.
This convention is stated in each result and does not use the deleted-fit
covariance in the Cook quadratic form.

Full random and residual coefficient covariance matrices and their changes
are retained in their original coefficient units. Unidentified covariance
coordinates should be assessed through implied covariance and recorded
boundary diagnostics, rather than interpreted as precisely estimated
individual parameters.
}
\subsection{Predictive studentization}{
The deleted fit predicts the held-out response using its full Gaussian
covariance, including random-effect and residual cross-covariance.
The indices \eqn{o} and \eqn{t} identify the remaining observed and deleted
cells, respectively. The cross-covariance is \eqn{K=V_{to}}, and
\eqn{C=(X_o'V_{oo}^{-1}X_o)^{-1}} is defined on an independent fixed basis.
The predictor, before restoring known offsets, is
\eqn{X_t\hat\beta+KV_{oo}^{-1}(y_o-X_o\hat\beta)}.
With \eqn{D=X_t-KV_{oo}^{-1}X_o}, its prediction-error covariance is
\eqn{V_{tt}-KV_{oo}^{-1}K'+DCD'}.

\code{rstudent} divides the observed-minus-predicted response by the
square root of the corresponding diagonal of this covariance. All observed
traits of a case are deleted together. The covariance parameters are fixed
at their deletion estimates for this calculation. No exact t reference
distribution is asserted for general mixed models. For ordinary Gaussian
regression fitted by REML this agrees with conventional externally
studentized residuals.
}
\subsection{Plots}{
The default panels show Cook's distance and the maximum absolute DFBETA
for each unit. \code{parm} selects individual signed DFBETA panels.
Unavailable and unrequested units retain their status and are counted in
the plot footer; no zero influence is substituted. If a panel has no finite
values it displays that fact explicitly. Plots annotate the largest
displayed values and invisibly return their exact points and annotations.
They do not alter the graphics layout or assign a scientific explanation
to an influential case or group.
}
}
\value{
\code{influence} returns class \code{st.lmm.influence}, containing an
all-unit \code{summary}, \code{inventory}, per-unit \code{results},
\code{requested} IDs, \code{fixed.change}, \code{dfbetas}, \code{Cook},
fitting controls, method and diagnostic definitions. The inventory maps
retained rows to stable unit IDs. Each result records row and response-cell
keys, status, warnings, fit diagnostics, lost strata, estimates, changes,
and a held-out prediction table and full prediction-error covariance when
available. Complete deletion fits are included only with
\code{keep.fits=TRUE}.

\code{summary} returns the all-unit table. \code{cooks.distance} returns
a named vector; \code{dfbeta} and \code{dfbetas} return matrices with one
row per inventory unit and one column per original fixed coefficient.
\code{rstudent} returns a matrix with one row per observation and one column
per trait, restores omitted
selected rows under \code{na.exclude}, and attaches the deletion-status
table. Originally missing responses and unavailable deletion diagnostics
remain \code{NA}.
}
\seealso{
\code{\link{st.lmm}}, \code{\link{st.lmm.diagnostics}},
\code{\link{st.lmm.refit}}, \code{\link{st.lmm.components}}
}
\keyword{models}
