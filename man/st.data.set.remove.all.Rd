% copyright (C) 2001 Frisch/Maurer
%
\name{st.data.set.remove.all}
\alias{st.data.set.remove.all}
\title{Remove All SelectionTools Marker Data Sets}
\description{
Releases all user-defined SelectionTools marker data sets and resets the built-in default data set.
}
\usage{
st.data.set.remove.all()
}
\details{
All native storage belonging to registered user-defined marker data sets is released and their registry entries are removed. The contents of \code{"default"} are reset; the built-in default name itself remains registered.
}
\value{
Invisibly returns the number of user-defined data sets removed. The built-in \code{"default"} data set remains registered internally but is empty. Because \code{st.data.set.list()} omits empty data sets, it returns a zero-row data frame if no other non-empty data sets remain.
}
