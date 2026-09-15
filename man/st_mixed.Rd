\name{st_mixed}
\alias{st_mixed}
\title{Fit a single-kernel mixed model using the SelectionTools native implementation}

\description{
Fit a single-kernel mixed model using the SelectionTools native implementation.
}

\usage{
st_mixed(y, id = NULL, Z, Z_ids = NULL, pool1 = NULL, pool2 = NULL,
         hybrid_id = NULL, include_sca = FALSE, ridge = 1e-8,
         method = c("std", "vanraden", "sommer"),
         coding = c("detect", "0 0.5 1", "-1 0 1", "0 1 2"), varcmp = NULL,
         maxcyc = 100L, tol = 1e-10, calcVC = 1,
         model_type = c("kernel", "marker"),
         bootstrap = list(nObs=NULL, nMarkers=NULL, nRep=1, seed=NULL))
}

\arguments{
   \item{y}{Response vector.}
   \item{id}{Identifier vector.}
   \item{Z}{Marker or design matrix.}
   \item{Z_ids}{Identifiers corresponding to rows of the marker matrix.}
   \item{pool1}{Function argument used by this operation.}
   \item{pool2}{Function argument used by this operation.}
   \item{hybrid_id}{Function argument used by this operation.}
   \item{include_sca}{Logical flag controlling inclusion of specific combining ability.}
   \item{ridge}{Ridge or numerical-stabilization parameter.}
   \item{method}{Method to use.}
   \item{coding}{Marker coding convention.}
   \item{varcmp}{Variance-component values.}
   \item{maxcyc}{Maximum number of fitting cycles.}
   \item{tol}{Convergence tolerance.}
   \item{calcVC}{Function argument used by this operation.}
   \item{model_type}{Function argument used by this operation.}
   \item{bootstrap}{Function argument used by this operation.}
}
\value{
A named \code{list} returned by the native mixed-model implementation, or \code{NULL} if the native routine cannot construct a result. The component set depends on the selected mode.

For the single-population mode the list contains \code{beta0} (intercept), \code{varcmp} (variance components), \code{ghat} (predicted genetic values), \code{yhat} (predicted response values), \code{missing_ids} (identifiers not matched to marker rows), \code{diag} (diagnostic counts and marker-coding information), and \code{retval} (compiled status code).

For hybrid modes the list contains \code{beta0}, \code{varcmp}, \code{yhat_obs} (predictions aligned to the observed-response vector), \code{diag}, \code{retval}, and \code{ridge}; a native implementation may additionally provide a \code{debug} component.
}
