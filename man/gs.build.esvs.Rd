\name{gs.build.esvs}
\alias{gs.build.esvs}
\title{Build estimation and validation subsets from a base data set}

\description{
Build estimation and validation subsets from a base data set.
}

\usage{
gs.build.esvs(es.size, base.set = "default", estimation.set = "none",
              validation.set = "none", es_m = NULL, es_p = NULL,
              vs_m = NULL, vs_p = NULL, auxfiles = FALSE)
}

\arguments{
   \item{es.size}{Size of the estimation set.}
   \item{base.set}{Name of the base data set.}
   \item{estimation.set}{Name of the estimation data set.}
   \item{validation.set}{Name of the validation data set.}
   \item{es_m}{Marker-data output file for the estimation set. Required when \code{auxfiles = TRUE}.}
   \item{es_p}{Phenotype-data output file for the estimation set. Required when \code{auxfiles = TRUE}.}
   \item{vs_m}{Marker-data output file for the validation set. Required when \code{auxfiles = TRUE}.}
   \item{vs_p}{Phenotype-data output file for the validation set. Required when \code{auxfiles = TRUE}.}
   \item{auxfiles}{Logical. If \code{TRUE}, the four persistent output files are written and \code{es_m}, \code{es_p}, \code{vs_m}, and \code{vs_p} must all be supplied.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of constructing the estimation and validation subsets. Validation failures may return \code{NULL}.
}
