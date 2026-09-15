\name{st.info}
\alias{st.info}
\title{Print or handle a SelectionTools information message}

\description{
Print or handle a SelectionTools information message.
}

\usage{
st.info(lev, msg)
}

\arguments{
   \item{lev}{Information or verbosity level.}
   \item{msg}{Message text.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of issuing the requested SelectionTools information message.
}
