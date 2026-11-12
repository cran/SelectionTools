\name{gs.cross.validation}
\alias{gs.cross.validation}
\title{gs.cross.validation()}
\description{
Splits a given data set into estimation and prediction set. Estimates marker effects
with the specified model using the estimation set and predicts
direct genomic values (DGV) for the built prediction set.
}
\usage{
gs.cross.validation(estimation.method, n.ts, n.runs, hsq = 0.8, alpha = 1,
                    maxiter = 100, precision = 0.0001,
                    out.filename.r = "crosval", out.filename.u = "effects",
                    auxfiles = TRUE, data.set = "default")
}
\arguments{
   \item{estimation.method}{character. Specifies the model: "rmlc", "rrwe", "rmlv", "rmlr", "rrwr" "bayes"}
   \item{n.ts}{Number of individuals in the training set}
   \item{n.runs}{Number of cross validation runs}
   \item{hsq}{Heritability of the analyzed trait. Used for "rrwe", "rrwr"}
   \item{alpha}{Function argument.}
   \item{maxiter}{Maximum number of iterations used for convergence of the REML
                estimator. Used if scheme "rmlc" or "rmlv"}
   \item{precision}{Gives the precision as second stop value for iteration next to maxiter. If
                  the difference of the values of lambda in two subsequent iterations
                  is less than the given precision, convergence is reached and the iteration will stop.
                  Used if scheme "rmlc" or "rmlv"}
   \item{out.filename.r}{Character string retained for compatibility. Cross-validation correlations are returned directly from memory.}
   \item{out.filename.u}{Character string retained for compatibility. Estimated marker effects are returned directly from memory.}
   \item{auxfiles}{Logical. If \code{TRUE}, both cross-validation result tables are returned directly from memory.}
   \item{data.set}{character. Name of the data set that is used by the function.}
}

%\details{
%}

\value{
Invisibly returns a \code{list} when \code{auxfiles = TRUE} and cross-validation succeeds. Component \code{cor} is a numeric matrix of cross-validation correlations and component \code{eff} is a numeric matrix of estimated effects. Otherwise \code{NULL}.
}

%\note{
%
%}
