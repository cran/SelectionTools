\name{gs.single.marker.aov}
\alias{gs.single.marker.aov}
\title{Performs single-marker analysis for a genomic-selection data set}

\description{
Performs single-marker analysis for a genomic-selection data set.
}

\usage{
gs.single.marker.aov(alpha = 1, out.filename = "sm.reg", auxfiles = TRUE,
                     data.set = "default")
}

\arguments{
   \item{alpha}{Method-specific tuning or significance parameter.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the single-marker ANOVA results when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
