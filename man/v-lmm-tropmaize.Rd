\name{v-lmm-tropmaize}
\alias{v-lmm-tropmaize}
\alias{tropmaize}
\docType{data}
\title{Constructed Tropical Maize Complete-Block Trial}
\description{
Synthetic plot observations for the 264 tropical maize lines in
\code{v-tropmaize-phe}, arranged in three complete blocks in one environment.
The field-only mixed-model genotype BLUEs reproduce the original phenotypic
values. The dataset illustrates field-trial and genomic analyses; it is
not the original plot-level experiment.
}
\usage{data("v-lmm-tropmaize", package = "SelectionTools")}
\format{A data frame named \code{tropmaize} with 792 rows and four columns:
\describe{
\item{\code{block}}{Complete-block factor with three levels, \code{1}--\code{3}. Each block contains all 264 genotypes.}
\item{\code{plot}}{Integer plot position, 1--264, within each block.}
\item{\code{gen}}{Genotype factor with 264 levels. Labels match the individual identifiers in \code{v-tropmaize-phe} and \code{v-tropmaize-vcf}.}
\item{\code{y}}{Positive numeric plot response on the original phenotype scale. Values are not rounded. No additional measurement unit is assigned.}
}
There are no missing observations. Each genotype occurs once in each block.
All plots belong to one environment, so no environment column is needed.
The \code{construction} attribute records the seed, design, nominal block
effects, residual standard deviation and centring rule.
The data file uses xz compression at level 9.
}
\details{
The construction uses seed 20261014. Nominal block effects are -0.10, 0 and
0.10. Independent normal plot errors with nominal standard deviation 0.18
are generated and centred within each genotype across the three blocks.
The centred errors and block effects are added to the original genotype
values. For genotypes whose deviations would be too large, all three
deviations are multiplied by the same factor to bound their absolute values
by 0.85 times the target mean. This retains the genotype mean and keeps all
responses positive. The bound and centring impose constraints on these
synthetic observations; they are illustrative construction choices, not
variance estimates or independent errors from the original study.
Plots order is independently randomized within each complete block.

The genotype arithmetic means across the three complete blocks reproduce
\code{v.tropmaize.phe$yy}, matched by \code{v.tropmaize.phe$ii}, to
floating-point precision. With fixed genotype effects, a random block
effect and a common residual variance, genotype BLUEs equal these arithmetic
means, regardless of the positive block and residual variance components.

Genomic predictions are not constrained to reproduce these genotype means.

The data represent one environment; no \code{env} or \code{rep} column
is stored.
}
\source{
Constructed for the SelectionTools linear mixed-model vignette from the
genotype identifiers and phenotypic values in \code{v-tropmaize-phe}.
The original phenotype and marker source is Crossa et al. (2010),
\emph{Genetics} 186, 713--724, \doi{10.1534/genetics.110.118521}.
}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lsm}},
\code{\link{gs.make.A}}, \code{\link{v-tropmaize-phe}},
\code{\link{v-tropmaize-vcf}}}
\keyword{datasets}
