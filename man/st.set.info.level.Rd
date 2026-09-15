\name{st.set.info.level}
\alias{st.set.info.level}
\title{Set info level used by SelectionTools}

\description{
Set info level used by SelectionTools.
}

\usage{
st.set.info.level(level)
}

\arguments{
   \item{level}{Function argument used by this operation.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of setting the SelectionTools information level.
}
