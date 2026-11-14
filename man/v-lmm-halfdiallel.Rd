\name{v-lmm-halfdiallel}
\alias{v-lmm-halfdiallel}
\alias{halfdiallel}
\docType{data}
\title{Maize Half Diallel for Shared Parental GCA and Cross SCA}
\description{
Sugar measurements for 21 unordered crosses among seven maize parents,
with two technical repetitions per cross. These are the small half-diallel
data used in the quantitative-genetics vignette of \pkg{lme4breeding}.
}
\usage{
data("v-lmm-halfdiallel", package = "SelectionTools")
}
\format{
Loads \code{halfdiallel}, a data frame with 42 rows and five variables:
\describe{
\item{\code{rep}}{Integer technical-repetition identifier, 1 or 2.}
\item{\code{male}}{Male-parent factor with all seven parent identifiers as levels.}
\item{\code{female}}{Female-parent factor with the same seven levels as male.}
\item{\code{cross}}{Factor identifying the 21 unordered parental crosses,
labelled as \code{P1:P2}, for example.}
\item{\code{sugar}}{Numeric sugar-content response in the source units.}
}
There are no missing values. All 21 unordered pairs of distinct parents
occur twice; self crosses and reciprocal observations are absent.
}
\details{
Response values and repetition labels are unchanged from \code{DT_halfdiallel}, distributed with
\pkg{enhancer} version 1.1.2 and used by \pkg{lme4breeding}. The object is
renamed for the SelectionTools vignette, with a source attribute added.
Both parental roles are factors with common levels P1 through P7. The
unordered cross factor is included, so the model needs no preprocessing.
The redundant source column \code{geno} is omitted: \code{cross} identifies
the same crosses and supplies the SCA term.
It is self-contained: neither source package is needed to load it.
}
\value{A data frame named \code{halfdiallel}, loaded into the selected
environment by \code{data}.}
\source{
Giovanny Covarrubias-Pazaran, \pkg{enhancer} 1.1.2,
\code{DT_halfdiallel}; quantitative-genetics vignette of
\pkg{lme4breeding}. The source package declares GPL version 2 or later.
}
\references{
Covarrubias-Pazaran G (2016) Genome-assisted prediction of quantitative
traits using the R package sommer. PLoS ONE 11:e0156744.
\doi{10.1371/journal.pone.0156744}
}
\keyword{datasets}
