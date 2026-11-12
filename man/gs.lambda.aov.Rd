\name{gs.lambda.aov}
\alias{gs.lambda.aov}
\title{Estimates or sets genomic-selection lambda values using the aov method}

\description{
Estimates or sets genomic-selection lambda values using the aov method.
}

\usage{
gs.lambda.aov(hsq, alpha = 1, out.filename = "lambda", auxfiles = TRUE,
              data.set = "default")
}

\arguments{
   \item{hsq}{Heritability parameter.}
   \item{alpha}{Method-specific tuning or significance parameter.}
   \item{out.filename}{Character string retained for compatibility. The result is returned directly from memory; this argument does not name an output file.}
   \item{auxfiles}{Logical. If \code{TRUE}, the result is returned directly from memory.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
Invisibly returns a numeric \code{matrix} containing the lambda results from the ANOVA method when \code{auxfiles = TRUE} and the compiled calculation succeeds; otherwise \code{NULL}.
}
