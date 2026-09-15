\name{mab.df}
\alias{mab.df}
\title{Helper used by the marker-assisted backcross simulation routines}

\description{
Resolves MAB file names against a configured directory while preserving
explicit absolute paths.
}

\usage{
mab.df(d, f)
}

\arguments{
   \item{d}{Directory used for a relative file name.}
   \item{f}{Non-empty file name or path. Absolute paths are returned unchanged.}
}
\value{
A character scalar. If \code{f} is absolute it is returned unchanged. If \code{f} is relative and \code{d} is non-empty, the two are combined with \code{file.path()}; otherwise \code{f} is returned unchanged.
}
