\name{st.lmm.model}
\alias{st.lmm.model}
\alias{model.frame.st.lmm}
\alias{terms.st.lmm}
\alias{weights.st.lmm}
\alias{update.st.lmm}
\title{Model Information and Updates for Gaussian Mixed Models}
\description{
Extracts retained model information or constructs a new fit from an
\code{st.lmm} object. The original fitted object is not modified.
}
\usage{
\method{model.frame}{st.lmm}(formula, ...)
\method{terms}{st.lmm}(x, ...)
\method{weights}{st.lmm}(object, ...)
\method{update}{st.lmm}(object, formula., ..., evaluate = TRUE)
}
\arguments{
\item{formula,x,object}{An \code{st.lmm} fitted object.}
\item{formula.}{An optional formula update, such as \code{. ~ . + x},
processed by \code{\link[stats]{update.formula}}. The named argument
\code{fixed} in \code{...} can instead replace the entire fixed formula.}
\item{evaluate}{A single nonmissing logical value. TRUE fits the updated
model; FALSE returns the call without evaluating additional arguments.}
\item{...}{For \code{update}, named arguments to \code{\link{st.lmm}}.
For the extraction methods, reserved for compatibility with the generics.}
}
\details{
The retained model frame contains the response, evaluated fixed-formula
variables, raw variables used by random and residual formulas, and supplied
weights and argument offsets. Fitted factor levels and contrast matrices are retained. The terms
object records the fitted prediction recipes for transformations such as
\code{poly} and \code{scale}. \code{model.response}, \code{model.offset}
and \code{model.weights} can extract their respective components from the
returned frame. Matrix responses remain matrices.

Under the default missing-response policy, a row with a missing response can
remain a prediction row. Its presence in \code{model.frame} does not imply
that it contributed to the likelihood. The fitted object's
\code{model$input.rows} and \code{model$row.names} retain row identities;
\code{nobs} reports the observed response count used by the fit. Under explicit
\code{na.omit} or \code{na.exclude}, the frame contains only retained rows;
the latter pads public response-level extraction without padding the frame or
engine matrices. The frame's \code{na.action} attribute identifies removed
positions in the selected row order.

The \code{model$input} record keeps the original data snapshot,
\code{selected.rows} (original row numbers, with repetitions allowed),
\code{selected.row.names}, and \code{retained.positions} within that selected
layout. \code{model$input.rows} gives the retained original row numbers.
\code{na.removed}, \code{weight.removed} and \code{common.complete.removed}
record distinct removal stages. \code{prediction.omitted} identifies retained
rows whose predictors were unusable under the default policy. These maps
refer to data rows; a multivariate likelihood may use several response cells
per retained row.

An update starts from the original input snapshot and fitted settings,
then applies its subset and missing-value policy again. It does not subset an
already selected frame. A changed formula can therefore recover rows whose
previously missing variables are no longer needed.
Changing or removing the caller's original data binding therefore does not
change a subsequent update. Supplying a new \code{data} argument explicitly
uses the new data. Data-derived subset, weight and offset expressions are
reevaluated there. An external training vector for subset, weights or offsets
requires an explicit replacement when changing data; it is not silently
recycled or matched only by length. Formula environments continue to resolve functions and
external objects referenced by the formulas; callers should keep such
objects available or include the required variables in the data frame.

Existing fitted factor coding is retained unless replaced through
\code{contrasts}. If an updated subset removes factor levels, the saved
contrast rows are matched to the remaining levels while retaining the fitted
contrast columns; resulting aliases are reported. Additional arguments are kept as expressions until the
returned call is evaluated, so data-dependent arguments such as
\code{weights = precision} are evaluated by \code{st.lmm} in its usual data
environment. The returned call uses the package-qualified fitting function
and can be evaluated after saving and reloading the model.

When the fitted object contains a subset search, an update changes the
specified full model. The search is repeated only when \code{selection} is
supplied explicitly; a warning explains this behavior. No selected candidate
automatically replaces the fitted model.
}
\value{
\code{model.frame} returns a data frame with a terms attribute;
\code{terms} returns a terms object; \code{weights} returns the retained
numeric precision vector, with zero elements for an unweighted fit.
\code{update} returns a new
\code{st.lmm} object, or an unevaluated call when \code{evaluate = FALSE}.
}
\seealso{\code{\link{st.lmm}}, \code{\link[stats]{model.frame}},
\code{\link[stats]{update}}}
\keyword{models}
