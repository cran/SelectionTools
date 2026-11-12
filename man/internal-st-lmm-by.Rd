\name{.st_lmm_by_prepare}
\alias{.st_lmm_by_prepare}
\alias{.st_lmm_by_task}
\alias{.st_lmm_by_collect}
\alias{.st_lmm_by_capture}
\title{Compiled Grouped Model Work Units}
\description{
Prepares independent grouped model fits, evaluates individual group tasks,
and collects their results. The condition callback captures errors and
warnings from compiled operations. These are implementation helpers for
\code{st.lmm.by}.
}
\usage{
.st_lmm_by_prepare(fixed, group, data, random = NULL, cov = ~units,
                   method = "REML", ddf = "Satterthwaite", type = 3L,
                   random.test = "none", control = list(), ucov = NULL,
                   weights = NULL, contrasts = NULL, subset = NULL,
                   na.action = NULL, offset = NULL, pool = FALSE)
.st_lmm_by_task(prepared, id)
.st_lmm_by_collect(prepared, results)
.st_lmm_by_capture(args)
}
\arguments{
\item{fixed, group, data, random, cov, method, ddf, type, random.test,
control, ucov, weights, contrasts, subset, na.action, offset, pool}{Model
specification, grouping, data selection and fitting arguments, as described
for \code{\link{st.lmm.by}}.}
\item{prepared}{The grouped task inventory returned by
\code{.st_lmm_by_prepare}.}
\item{id}{One integer group ID from the prepared inventory.}
\item{results}{A list containing exactly one completed task record for
each group ID in the prepared inventory.}
\item{args}{An internal compiled method request.}
}
\value{
The preparation helper returns an \code{st.lmm.by.prepared} inventory
containing group IDs, row mappings and task specifications. The task helper
returns one group result with fitting status and captured conditions.
Collection returns an \code{st.lmm.by} object. The condition callback
returns \code{value}, \code{error} and \code{warnings}.
}
\keyword{internal}
