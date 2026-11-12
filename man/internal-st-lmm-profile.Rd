\name{.st_lmm_profile_prepare}
\alias{.st_lmm_profile_prepare}
\alias{.st_lmm_profile_task}
\title{Compiled Likelihood Profile Work Units}
\description{
Prepares a common fitted-model baseline and evaluates an individual
likelihood profile parameter. These are implementation helpers for
\code{st.lmm.profile}.
}
\usage{
.st_lmm_profile_prepare(object, level = 0.95, control = list())
.st_lmm_profile_task(prepared, parm)
}
\arguments{
\item{object}{A fitted \code{st.lmm} model.}
\item{level}{The confidence level, as described for
\code{\link{st.lmm.profile}}.}
\item{control}{Named profiling and fitting controls, as described for
\code{\link{st.lmm.profile}}.}
\item{prepared}{The baseline and parameter inventory returned by
\code{.st_lmm_profile_prepare}.}
\item{parm}{One parameter selected by name or index from the prepared
inventory.}
}
\value{
The preparation helper returns the baseline model, parameter inventory,
confidence level and controls. The task helper returns the likelihood
profile result for the selected parameter.
}
\keyword{internal}
