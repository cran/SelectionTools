\name{info.cat}
\alias{info.cat}
\title{Information-level-controlled status output}

\description{
A status message is signaled when the requested information level is enabled.
The message is signaled with \code{message()} and can therefore also be
suppressed with \code{suppressMessages()}.
}

\usage{
info.cat(lev, msg)
}

\arguments{
   \item{lev}{Information or verbosity level.}
   \item{msg}{Message text.}
}
\value{
The function returns \code{NULL} invisibly. Message text is evaluated only when the requested information level is enabled.
}
