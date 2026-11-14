\name{st.set.openblas.threads}
\alias{st.set.openblas.threads}
\title{Sets the linked BLAS thread count}

\description{
Sets the linked BLAS thread count.
}

\usage{
st.set.openblas.threads(active)
}

\arguments{
   \item{active}{One finite whole number between 1 and 256.}
}
\details{
The package and OpenMP thread budgets remain unchanged. The linked BLAS
backend must expose supported thread controls and permit an independent
setting. An unavailable control or incompatible shared OpenMP setting
produces an error.
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of setting the OpenBLAS thread count.
}
