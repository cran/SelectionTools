\name{st.recode.ref.2}
\alias{st.recode.ref.2}
\title{Recode reference-marker information using the package-specific convention}

\description{
Recode reference-marker information using the package-specific convention.
}

\usage{
st.recode.ref.2(r1 = 0, r2 = 1, descending = FALSE,
                data.set = "default")
}

\arguments{
   \item{r1}{Function argument used by this operation.}
   \item{r2}{Function argument used by this operation.}
   \item{descending}{Function argument used by this operation.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of recoding reference-marker information for two parents. Validation failures may return \code{NULL}.
}
