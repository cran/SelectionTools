\name{st.start.timer}
\alias{st.start.timer}
\title{Starts the SelectionTools timing helper}

\description{
Starts the SelectionTools timing helper.
}

\usage{
st.start.timer(depth=1)
}

\arguments{
   \item{depth}{Timer depth, a finite whole number from 1 through 5.}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of starting the SelectionTools timing helper.
}
