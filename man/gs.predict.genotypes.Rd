\name{gs.predict.genotypes}
\alias{gs.predict.genotypes}
\title{gs.predict.genotypes()}
\description{
Predicts direct genomic values (DGV) for a given set of individuals in the prediction set.
Works only after the marker effects were estimated.
}
\usage{
gs.predict.genotypes(training.set ="default", prediction.set ="default",
                     out.filename = "breeding.values", auxfiles = TRUE)
}
\arguments{
   \item{training.set}{character. Name of the data set that was used for estimating
                        marker effects}
   \item{prediction.set}{character. Name of the data set that is used for predicting DGVs}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
}

%\details{
%
%}

\value{
Invisibly returns a \code{data.frame} with prediction results when \code{auxfiles = TRUE} and prediction succeeds. It contains individual identifiers and predicted genomic values, and also observed phenotypes when those are present in the generated result. Otherwise \code{NULL}.
}

%\note{
%
%}
