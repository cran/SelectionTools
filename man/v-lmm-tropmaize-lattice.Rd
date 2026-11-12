\name{v-lmm-tropmaize-lattice}
\alias{v-lmm-tropmaize-lattice}
\alias{tropmaize.lattice}
\docType{data}
\title{Tropical Maize Alpha-Lattice Trials in Five Environments}
\description{
A constructed tropical maize dataset with 264 genotypes evaluated in
three replications at each of five environments. Each replication is an
alpha-lattice with 22 incomplete blocks of 12 plots.
}
\usage{data("v-lmm-tropmaize-lattice", package = "SelectionTools")}
\format{
A data frame named \code{tropmaize.lattice} with 3960 rows and six columns:
\describe{
\item{\code{env}}{Environment factor with levels \code{E1} to \code{E5}.}
\item{\code{rep}}{Replication factor with three levels within each environment.}
\item{\code{block}}{Incomplete-block factor with 22 levels within each replication and environment.}
\item{\code{plot}}{Integer plot position, 1 to 12, within each incomplete block.}
\item{\code{gen}}{Genotype factor with 264 levels. Identifiers match those in \code{v-tropmaize-phe} and \code{v-tropmaize-vcf}.}
\item{\code{y}}{Positive numeric plot response on the scale of the supplied maize phenotypes. No additional measurement unit is assigned.}
}
There are no missing observations. Every genotype occurs once in each
replication, giving 15 observations per genotype and 792 per environment.
The \code{construction} attribute retains the seed and construction settings.
}
\details{
Replication and block labels are reused across environments. A replication
is identified by \code{env:rep}, and an incomplete block by
\code{env:rep:block}. Each pair of distinct genotypes occurs together in
an incomplete block at most once within an environment. The treatment
information matrix has rank 263 in every environment.

This is the original five-environment lattice dataset. It is available
separately from \code{v-lmm-tropmaize}, the three-complete-block dataset
in one environment. Loading this topic creates \code{tropmaize.lattice}.

The responses are synthetic, not observations from the original field
experiment. The construction uses common genotype-to-block incidence across
environments and centred contributions so that field-model genotype BLUEs
recover the supplied phenotype values. These imposed constraints should be
considered when interpreting analyses or using the data for benchmarks.
The original observations and construction attributes are retained unchanged.

}
\source{
Constructed for SelectionTools mixed-model examples from the genotype
identifiers and phenotypic values in \code{v-tropmaize-phe}.
}
\seealso{\code{\link{v-lmm-tropmaize}}, \code{\link{st.lmm}},
\code{\link{st.lsm}}, \code{\link{gs.make.A}},
\code{\link{v-tropmaize-phe}}, \code{\link{v-tropmaize-vcf}}}

\keyword{datasets}
