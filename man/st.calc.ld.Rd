\name{st.calc.ld}
\alias{st.calc.ld}
\title{Calculates Linkage Disequilibrium from Marker Data}

\description{
Calculates pairwise linkage disequilibrium within chromosomes. LD can either
be calculated directly from the two stored homologues or, for biallelic
markers, after estimation of unknown two-locus gametic phase.
}

\usage{
st.calc.ld(ld.measure = "r2", auxfiles = TRUE, data.set = "default")
}

\arguments{
  \item{ld.measure}{Character string selecting the LD measure. Supported
  values are \code{"r2"}, \code{"Dp"}, \code{"r2-estimate-phases"},
  \code{"Dp-estimate-phases"}, and \code{"r-estimate-phases"}. The last
  option returns the signed gametic correlation after phase estimation.}
  \item{auxfiles}{Logical. If \code{TRUE}, the pairwise results are returned directly
  from memory as a data frame. If \code{FALSE}, LD is calculated and retained internally in
  the data set but the pairwise table is not transferred to R.}
  \item{data.set}{Character string naming the SelectionTools marker data set.}
}

\details{
LD is calculated for every marker pair on the same chromosome.
Chromosome groups and marker order are obtained from the data set's map.

With \code{ld.measure="r2"} and \code{ld.measure="Dp"}, the marker
representation is treated as phased. The two allele arrays represent the two
stored homologues, and alleles at two loci from the same homologue are counted
together as haplotypes. Thus, for heterozygous material the stored allele order
must represent the homologues correctly. For completely homozygous lines the
phase is immaterial. These two direct calculations support multiallelic loci.

For allele \eqn{x} at the first locus and allele \eqn{y} at the second locus,
\eqn{f_x} and \eqn{f_y} denote the allele frequencies and \eqn{f_{xy}} denotes
the haplotype frequency. The gametic disequilibrium coefficient is
\deqn{D_{xy}=f_{xy}-f_x f_y.}
For \code{ld.measure="r2"}, the allele-pair value is
\deqn{r^2_{xy}=\frac{D_{xy}^2}{f_x(1-f_x)f_y(1-f_y)}}
and the multiallelic value returned by SelectionTools is
\deqn{r^2=\sum_x\sum_y f_x f_y r^2_{xy}.}

For \code{ld.measure="Dp"}, the allele-pair value is the absolute normalized
disequilibrium
\deqn{D'_{xy}=\left|\frac{D_{xy}}{D_{max}}\right|,}
where
\deqn{D_{max}=\min(f_x f_y,(1-f_x)(1-f_y)) \quad \mbox{if } D_{xy}<0}
and
\deqn{D_{max}=\min(f_x(1-f_y),(1-f_x)f_y) \quad \mbox{if } D_{xy}\ge 0.}
The multiallelic value is
\deqn{D'=\sum_x\sum_y f_x f_y D'_{xy}.}

The options \code{"r2-estimate-phases"},
\code{"Dp-estimate-phases"}, and \code{"r-estimate-phases"} are intended
for heterozygous material with
unknown gametic phase. They require exactly two biological alleles at every
marker in the selected data set. The missing-allele state is not a biological
allele. Individuals with missing data at either locus of a marker pair are
omitted from that pair's two-locus genotype table.

For the phase-estimating modes, the nine observable genotype classes at two
biallelic loci are counted without using the order of the two stored homologues
as phase information. All genotype classes except the double heterozygote have
an unambiguous pair of haplotypes. The quantities \eqn{h_{00}}, \eqn{h_{01}},
\eqn{h_{10}}, and \eqn{h_{11}} denote the current haplotype-frequency
estimates. In the expectation step, the probability that a double heterozygote
has the coupling diplotype is
\deqn{q=\frac{h_{00}h_{11}}
             {h_{00}h_{11}+h_{01}h_{10}}.}
The ambiguous individuals are split between the coupling and repulsion
diplotypes according to \eqn{q}; the four haplotype frequencies are then
re-estimated from the resulting expected haplotype counts. These expectation
and maximization steps are repeated to convergence. SelectionTools uses several
deterministic starting values and retains the converged solution with the
greatest observed-data likelihood. The corresponding \eqn{r}, \eqn{r^2}, or \eqn{D'}
is calculated from the estimated haplotype frequencies. The estimates maximize
the observed-data likelihood for the two-locus phase-ambiguity problem;
Clayton and Leung (2007) discuss haplotype estimation for association data.

For \code{"r-estimate-phases"}, the signed statistic is
\deqn{r=\frac{h_{00}-p_0q_0}{\sqrt{p_0(1-p_0)q_0(1-q_0)}}.}
Here allele zero is the first biological allele in the marker's sorted
internal allele-code order, and \eqn{p_0} and \eqn{q_0} are its two marginal
frequencies among complete pairs. Positive values indicate coupling of these
two alleles; exchanging the allele labels at one locus reverses the sign.
The statistic is a gametic correlation estimated from the two-locus
likelihood, not a correlation of observed diploid allele dosages. Reversing
the order of the two alleles within an individual's unphased genotype has no
effect.

Signed \eqn{r} is missing if no complete genotype pairs remain, a locus
becomes monomorphic among complete pairs, the EM calculation does not yield
a converged solution, or equally supported likelihood maxima have opposite
signs. The last case can occur when all complete observations are double
heterozygotes; their phase does not identify a sign. The existing \eqn{r^2}
and \eqn{D'} phase-estimation paths retain their zero value for a pair with
no usable variation or no converged solution.

For the direct \code{"r2"} and \code{"Dp"} calculations, the missing allele
is excluded from the allele-pair sums, while the frequency denominator is
twice the number of individuals in the data set. Marker filtering for
missingness and polymorphism can therefore be useful before LD calculation.

When \code{auxfiles=FALSE}, the internally stored LD values can be used by
\code{st.def.hblocks} without constructing a potentially very large R data
frame. This can reduce result-transfer time for large data sets.

A general review of LD measures and their dependence on allele frequencies
is given by Hedrick (1987), \emph{Genetics} 117:331--341.
}

\value{
If \code{auxfiles=TRUE} and calculation succeeds, invisibly returns a
\code{data.frame} with columns \code{Chrom}, \code{Locus1}, \code{Locus2},
\code{Name1}, \code{Name2}, and \code{LD}. \code{Chrom} is the one-based
chromosome index in the selected data set. \code{Locus1} and \code{Locus2}
are one-based marker positions within that chromosome. If \code{auxfiles=FALSE},
invisibly returns \code{NULL}; the calculated LD values remain stored in the
selected data set for subsequent SelectionTools operations. If a
phase-estimating measure is requested and the data set is not biallelic, the
calculation is rejected.
}

\references{
Clayton D, Leung HT (2007) An R package for analysis of whole-genome
association studies. Human Heredity 64:45--51.

Hedrick PW (1987) Gametic disequilibrium measures: proceed with caution.
Genetics 117:331--341.
}

