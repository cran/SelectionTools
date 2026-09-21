\name{st.calc.ld}
\alias{st.calc.ld}
\title{Calculate Linkage Disequilibrium from Marker Data}

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
  values are \code{"r2"}, \code{"Dp"}, \code{"r2-estimate-phases"}, and
  \code{"Dp-estimate-phases"}.}
  \item{auxfiles}{Logical. If \code{TRUE}, the pairwise results are written to
  an R-session temporary file, read back by the R wrapper, and returned as a
  data frame. If \code{FALSE}, LD is calculated and retained internally in
  the data set but the pairwise table is not transferred to R.}
  \item{data.set}{Character string naming the SelectionTools marker data set.}
}

\details{
LD is calculated for every marker pair on the same chromosome.

With \code{ld.measure="r2"} and \code{ld.measure="Dp"}, the marker
representation is treated as phased. The two allele arrays represent the two
stored homologues, and alleles at two loci from the same homologue are counted
together as haplotypes. Thus, for heterozygous material the stored allele order
must represent the homologues correctly. For completely homozygous lines the
phase is immaterial. These two direct calculations support multiallelic loci.

For allele \eqn{x} at the first locus and allele \eqn{y} at the second locus,
let \eqn{f_x} and \eqn{f_y} be the allele frequencies and \eqn{f_{xy}} the
haplotype frequency. The gametic disequilibrium coefficient is
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

The options \code{"r2-estimate-phases"} and
\code{"Dp-estimate-phases"} are intended for heterozygous material with
unknown gametic phase. They require exactly two biological alleles at every
marker in the selected data set. The missing-allele state is not a biological
allele. Individuals with missing data at either locus of a marker pair are
omitted from that pair's two-locus genotype table.

For the phase-estimating modes, the nine observable genotype classes at two
biallelic loci are counted without using the order of the two stored homologues
as phase information. All genotype classes except the double heterozygote have
an unambiguous pair of haplotypes. Let \eqn{h_{00}}, \eqn{h_{01}},
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
greatest observed-data likelihood. The corresponding \eqn{r^2} or \eqn{D'}
is calculated from the estimated haplotype frequencies. This is a native
maximum-likelihood implementation for the two-locus phase-ambiguity problem;
Clayton and Leung (2007) discuss haplotype estimation for association data.

For the direct \code{"r2"} and \code{"Dp"} calculations, the missing allele
is excluded from the allele-pair sums, while the frequency denominator is
twice the number of individuals in the data set. Marker filtering for
missingness and polymorphism can therefore be useful before LD calculation.

When \code{auxfiles=FALSE}, the internally stored LD values can be used by
\code{st.def.hblocks} without constructing a potentially very large R data
frame. This can reduce transfer and file-reading time for large data sets.

For a general review of LD measures and their dependence on allele frequencies,
see Hedrick (1987), \emph{Genetics} 117:331--341.
}

\value{
If \code{auxfiles=TRUE} and calculation succeeds, invisibly returns a
\code{data.frame} with columns \code{Chrom}, \code{Locus1}, \code{Locus2},
\code{Name1}, \code{Name2}, and \code{LD}. \code{Locus1} and \code{Locus2}
are one-based marker positions within the chromosome. If \code{auxfiles=FALSE},
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

\seealso{
\code{st.LDplot.ld}, \code{st.LDplot.map}, \code{st.LDheatmap},
\code{st.def.hblocks}
}
