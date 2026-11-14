\name{st.stop.timer}
\alias{st.stop.timer}
\title{Stops the SelectionTools timing helper}

\description{
Stops the SelectionTools timing helper.
}

\usage{
st.stop.timer(depth=1, info.level=0)
}

\arguments{
   \item{depth}{Timer depth, a finite whole number from 1 through 5.}
   \item{info.level}{Function argument used by this operation.}
}
\value{
An invisible implementation-level \code{list} returned by the final information-message call after elapsed time is calculated and reported.
}
