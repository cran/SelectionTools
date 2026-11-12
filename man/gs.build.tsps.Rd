\name{gs.build.tsps}
\alias{gs.build.tsps}
\title{Builds training and prediction subsets from a base data set}

\description{
Builds training and prediction subsets from a base data set.
}

\usage{
gs.build.tsps(ts.size, base.set = "default", training.set = "none",
              prediction.set = "none", es_m = NULL, es_p = NULL,
              vs_m = NULL, vs_p = NULL, auxfiles = FALSE)
}

\arguments{
   \item{ts.size}{Size of the training set.}
   \item{base.set}{Name of the base data set.}
   \item{training.set}{Name of the training data set.}
   \item{prediction.set}{Name of the prediction data set.}
   \item{es_m}{Marker-data output file for the training set. Required when \code{auxfiles = TRUE}.}
   \item{es_p}{Phenotype-data output file for the training set. Required when \code{auxfiles = TRUE}.}
   \item{vs_m}{Marker-data output file for the prediction set. Required when \code{auxfiles = TRUE}.}
   \item{vs_p}{Phenotype-data output file for the prediction set. Required when \code{auxfiles = TRUE}.}
   \item{auxfiles}{Logical. If \code{TRUE}, the four persistent output files are written and \code{es_m}, \code{es_p}, \code{vs_m}, and \code{vs_p} must all be supplied.}
}
\details{
The two output sets must be distinct from each other and from the source.
Both must contain at least one individual. Original allele codes, homologue
order, phase state, linkage map and the selected phenotypic values are retained.
Each output receives current marker statistics and automatically selected
4-, 8- or 16-bit allele storage. Both outputs are prepared and checked before
either existing output is replaced.
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of constructing the training and prediction subsets. Validation failures may return \code{NULL}.
}
