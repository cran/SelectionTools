\name{gs.w.aov}
\alias{gs.w.aov}
\title{Calculates weights or statistics used by the genomic-selection ANOVA method}

\description{
Calculates weights or statistics used by the genomic-selection ANOVA method.
}

\usage{
gs.w.aov(alpha = 1, out.filename = "W", auxfiles = TRUE, data.set = "default")
}

\arguments{
   \item{alpha}{Method-specific tuning or significance parameter.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the ANOVA weights/statistics when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
