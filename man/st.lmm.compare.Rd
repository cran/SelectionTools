\name{st.lmm.compare}
\alias{st.lmm.compare}
\alias{drop1.st.lmm}
\title{Likelihood Comparisons and Fixed-Term Deletion for Mixed Models}
\description{
Compares fitted univariate Gaussian mixed models using likelihoods and
information criteria, or refits models after deleting eligible fixed terms.
Likelihood-ratio inference is supplied only when nesting has been established.
}
\usage{
\method{drop1}{st.lmm}(object, scope, test = "none", k = 2,
                      trace = FALSE, ..., refit = TRUE)
}
\arguments{
\item{object}{A converged univariate \code{st.lmm} model fitted by ML or REML.}
\item{scope}{Optional character vector of fixed-term names or a one-sided
formula identifying those terms. By default, all terms eligible under the
marginality hierarchy are considered. An interaction prevents deletion of
its constituent lower-order terms. The intercept is not a deletion term.}
\item{test}{\code{"none"} reports information criteria; \code{"Chisq"}
also reports likelihood-ratio statistics and approximate chi-squared
probabilities. The multi-model \code{anova} interface defaults to
\code{"Chisq"}; \code{drop1} defaults to \code{"none"}.}
\item{k}{Finite nonnegative information-criterion penalty per estimated
parameter. The default 2 gives AIC; \code{log(nobs(object))} gives BIC.}
\item{trace}{Logical; TRUE prints the term being deleted.}
\item{refit}{Logical. For multi-model \code{anova}, TRUE refits every REML
model by ML. FALSE keeps the fitted method and permits REML comparison only
with identical fixed design matrices. Fixed-term deletion requires ML, so a
REML object is refitted when TRUE and rejected when FALSE.}
\item{...}{Reserved for \code{drop1}. For multi-model \code{anova},
additional fitted models are unnamed arguments.}
}
\details{
\code{anova(m1, m2, ...)} compares ML or REML \code{st.lmm} objects.
The \code{effects="fixed"} and \code{effects="random"} selections return
the corresponding tests cached during fitting. Both \code{effects="both"}
and \code{effects="all"} return the combined table. Moments-model comparisons retain
their separate exact-regression behavior.

All compared models must retain identical response values, observation
identities and order, precision weights, and offsets. This includes the
layout of any retained missing-response prediction rows. Models fitted to
different observations cannot be compared through this interface. Every
supplied fit must have a finite likelihood and report convergence.

The default ML refitting makes comparisons of different fixed effects
meaningful. Refitting is reported in the table note and its
\code{ML.refitted} attribute. The supplied objects are not changed.
With \code{refit=FALSE}, ML and REML cannot be mixed. REML comparisons also
require the same numeric fixed design matrix: using different contrast
bases can change the additive constant in a REML likelihood even
when the fixed model spaces coincide.

Models are ordered by their number of free parameters. This count is the
estimable fixed-effect rank plus the nominal number of free covariance
parameters. A fitted variance of zero does not remove its parameter from
this count. AIC and BIC use this count and the observed response count.
The \code{input.order} attribute maps table rows to the supplied order.
The optional \code{model.names} argument to \code{anova} supplies unique
nonempty labels; otherwise concise, distinct model expressions are used
with numbered labels as a fallback.

For consecutive models with a positive parameter-count difference, nesting
is checked using the fixed design spaces and the random covariance
families. Supported random comparisons include removal of free random
blocks, scalar random coefficients within a larger unstructured block,
and independent coefficients within a corresponding unstructured block.
The homogeneous-diagonal, compound-symmetry and AR(1) families are compared
only when their covariance constraints establish inclusion. The residual
covariance families and physical bases must agree. For unproved nesting,
the table retains information criteria and reports unavailable LRT
inference with an explanation. Equal parameter counts also give no LRT
probability. No test is inferred from parameter counts alone.

Covariance-parameter LRTs use an approximate chi-squared distribution.
Variance parameters on the boundary of their allowable range can require
a different calibration; such comparisons have an explicit table note.
For supported whole-random-term tests, \code{st.lmm(random.test="bootstrap")}
provides the existing parametric-bootstrap alternative. A fixed-term
likelihood-ratio test differs from a cached Wald F test and does not use
Satterthwaite or Kenward--Roger denominator degrees of freedom.

\code{drop1} removes one complete eligible formula term per fit. A factor
term can therefore remove several estimable coefficients. The same
observations, offsets, weights and covariance specification must remain
in use. Dropping an aliased term can leave the estimable rank unchanged;
the table then gives no invented degrees of freedom or probability.
An unsuccessful reduced fit is retained in the table with an explanatory
\code{status} and unavailable inference. An incompatible observation
layout is an error. Models are reconstructed from retained input and
formula metadata rather than a mutable caller data binding.
}
\value{
A data frame using the package's table printer. Multi-model \code{anova}
reports \code{model}, \code{npar}, \code{AIC}, \code{BIC},
\code{logLik}, \code{deviance}, \code{Chisq}, \code{DF},
\code{p} and \code{comparison}. The \code{comparison} column
identifies the reference row, nested tests, or unavailable inference.

\code{drop1} reports \code{term}, \code{DF}, \code{AIC},
\code{LRT}, \code{p} and \code{status}; the first row is
\code{<none>}, the full model. Its \code{AIC} column uses the supplied
penalty \code{k}. The printer labels the column AIC for \code{k=2},
BIC for \code{k=log(nobs(object))}, and Criterion for other penalties.
The \code{method} and \code{ML.refitted} attributes
record the estimation convention.
}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lmm.model}},
\code{\link[stats]{anova}}, \code{\link[stats]{drop1}}}
\keyword{models}
