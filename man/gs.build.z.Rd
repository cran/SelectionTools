\name{gs.build.Z}
\alias{gs.build.Z}
\title{gs.build.Z()}
\description{
  The function \emph{gs.build.Z} builds the design (Z) matrix out of the marker data
  for genomic selection.
}
\usage{
gs.build.Z(out.filename = "Z.matrix", auxfiles = FALSE, data.set = "default")
}
\arguments{
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{character. Name of the data set that is used by the function}
}
\value{
Invisibly returns a \code{matrix} containing the generated marker/design matrix \code{Z} when \code{auxfiles = TRUE} and generation succeeds; otherwise \code{NULL}.
}

%\details{}

