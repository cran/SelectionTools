\name{info.cat}
\alias{info.cat}
\title{Emit information-level-controlled status output}

\description{
Emit a status message when the requested information level is enabled.
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
\code{NULL}. The function is called for its side effect of conditionally signaling a status message according to the current information level.
}
