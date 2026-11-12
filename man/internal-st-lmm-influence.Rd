\name{.st_lmm_influence_inventory}
\alias{.st_lmm_influence_inventory}
\alias{.st_lmm_delete_unit}
\alias{.st_lmm_influence_collect}
\alias{.st_lmm_diagnostic_capture}
\title{Compiled Deletion Diagnostic Work Units}
\description{
Constructs an exact fitted-model deletion inventory, evaluates one unit,
and collects unit results without discarding failures. The condition
callback captures errors and warnings from compiled fitting and diagnostic
operations. These are implementation helpers rather than a separate model
fitting interface.
}
\usage{
.st_lmm_influence_inventory(object, group = NULL)
.st_lmm_delete_unit(object, inventory, id, control = list(),
                    keep.fit = FALSE)
.st_lmm_influence_collect(object, inventory, results, control = list())
.st_lmm_diagnostic_capture(object, action, request)
}
\arguments{
\item{object}{A fitted \code{st.lmm} model.}
\item{group}{Deletion groups as described for \code{influence.st.lmm}.}
\item{inventory}{An inventory returned by the matching helper. It contains
all unit IDs and an exact serializable numeric/structural model identity.}
\item{id}{One integer inventory unit ID.}
\item{control}{Named fitting-control overrides.}
\item{keep.fit}{Retain the complete deletion model when TRUE.}
\item{results}{A list in inventory-ID order, with one slot per unit.
NULL slots are retained as not requested. Present records must belong to
the fitted model and their stated unit.}
\item{action}{An internal fitting or diagnostic operation.}
\item{request}{Its compiled request object.}
}
\value{
The inventory helper returns unit keys and row membership. The per-unit
helper returns a deletion record, including status and conditions.
Collection returns an \code{st.lmm.influence} object. The condition callback
returns \code{value}, \code{error} and \code{warnings}.
}
\seealso{\code{\link{influence.st.lmm}}}
\keyword{internal}
