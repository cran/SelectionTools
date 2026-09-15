\name{st.mxd}
\alias{st.mxd}
\title{Fit the package multi-kernel mixed-model routine}

\description{
Fit the package multi-kernel mixed-model routine.
}

\usage{
st.mxd(y, X, Z = NULL, V = NULL, estimates = FALSE, maxcyc = 100,
       precision = 1e-4)
}

\arguments{
   \item{y}{Response vector.}
   \item{X}{Fixed-effect design matrix.}
   \item{Z}{Marker or design matrix.}
   \item{V}{Function argument used by this operation.}
   \item{estimates}{Function argument used by this operation.}
   \item{maxcyc}{Maximum number of fitting cycles.}
   \item{precision}{Convergence tolerance.}
}
\value{
A named \code{list} on successful fitting, with components:

  \item{sigma}{Numeric vector of variance-component estimates.}
  \item{b}{Fixed-effect estimates when \code{estimates = TRUE}; otherwise the scalar initialization value used by the routine.}
  \item{u}{Random-effect estimates for \code{Z}-based models when \code{estimates = TRUE}; otherwise the scalar initialization value used by the routine.}

Returns \code{NULL} when input validation or the compiled fit reports failure.
}
