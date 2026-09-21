\name{show.phase}
\alias{show.phase}
\title{Display a Compact Statistical Phasing Summary}
\description{
Displays the main diagnostics from \code{st.phase()} together with switch-error
statistics from \code{st.switch.error()}.
}
\usage{
show.phase(label, result, score)
}
\arguments{
  \item{label}{Character string printed above the summary.}
  \item{result}{An object returned by \code{st.phase()}.}
  \item{score}{The one-column data frame returned by \code{st.switch.error()}.}
}
\value{
Returns \code{result} invisibly. The function is called for the side effect of
printing a compact phasing summary.
}
\seealso{
\code{st.phase}, \code{st.switch.error}, \code{st.destroy.phase}
}
