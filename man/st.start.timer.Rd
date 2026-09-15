\name{st.start.timer}
\alias{st.start.timer}
\title{Start or stop the SelectionTools timing helper}

\description{
Start or stop the SelectionTools timing helper.
}

\usage{
st.start.timer(depth=1)
}

\arguments{
   \item{depth}{Function argument used by this operation.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of starting the SelectionTools timing helper.
}
