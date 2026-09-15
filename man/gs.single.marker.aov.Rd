\name{gs.single.marker.aov}
\alias{gs.single.marker.aov}
\title{Perform single-marker analysis for a genomic-selection data set}

\description{
Perform single-marker analysis for a genomic-selection data set.
}

\usage{
gs.single.marker.aov(alpha = 1, out.filename = "sm.reg", auxfiles = TRUE,
                     data.set = "default")
}

\arguments{
   \item{alpha}{Method-specific tuning or significance parameter.}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the single-marker ANOVA results when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
