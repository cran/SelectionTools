\name{st.set.openblas.threads}
\alias{st.set.openblas.threads}
\title{Set openblas threads used by SelectionTools}

\description{
Set openblas threads used by SelectionTools.
}

\usage{
st.set.openblas.threads(active)
}

\arguments{
   \item{active}{Function argument used by this operation.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of setting the OpenBLAS thread count.
}
