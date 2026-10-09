\name{st.phase}
\alias{st.phase}
\title{Cohort and Reference-Assisted Phasing of Diploid Genotypes}
\description{
Infers haplotypes at existing markers from the cohort or from fixed reference
haplotypes. Missing genotypes at typed markers can optionally be filled.
Observed unordered genotypes are preserved. Automatic window selection uses
the genetic map and the spacing of common markers.
}
\usage{
st.phase(data.set = "default", ref = NULL, phasing.options = list(),
         impute = FALSE, commit = TRUE)
}
\arguments{
  \item{data.set}{Name of an existing SelectionTools marker data set with a
  linkage map. Genetic positions must be in centiMorgans. All individuals and
  all chromosomes in this data set are phased.}
  \item{ref}{Reference panel, either the name of a loaded SelectionTools
  data set or an \code{STvcf} object. With \code{NULL}, the cohort supplies
  donor haplotypes. Reference haplotypes remain fixed. See the reference
  requirements below. Default: \code{NULL}.}
  \item{phasing.options}{Uniquely named list of controls described below.
  Omitted entries use their defaults. Hyphens in option names are accepted as
  alternatives to dots, for example \code{"phase-states"}.
  Unknown or duplicate entries are errors.}
  \item{impute}{Logical scalar. If \code{FALSE}, preserve missing target
  genotypes. If \code{TRUE}, retain inferred hard calls at missing genotypes
  of existing target markers. Observed unordered genotypes are preserved
  under both policies. Default: \code{FALSE}.}
  \item{commit}{Logical scalar. If \code{TRUE}, write the completed result to
  \code{data.set} and mark it phased. If \code{FALSE}, return the result
  without changing the data set. Default: \code{TRUE}.}
}
\details{
\strong{Cohort and reference modes.}
With \code{ref = NULL}, the target cohort supplies the donor
haplotypes. With a reference panel, both target and reference haplotypes can
serve as donors, but only target individuals are updated. References remain
fixed during initialization, common-marker iterations, rare-marker
reconstruction and window processing. EM estimation uses targets; target
relatedness exclusions do not remove reference donors. This is a general
copying model: there is no F2-specific ancestry model or pedigree option.

The implementation selects common markers, initializes haplotypes, searches
for donors with the positional Burrows--Wheeler transform, iterates a reduced
haplotype-copying hidden Markov model, and reconstructs remaining rare-marker
alleles. Chromosomes and overlapping windows are processed in genetic-map
order. Sample work within several stages runs in parallel. More threads do
not remove the sequential dependencies between stages.

\strong{Input requirements.}
There must be 2--5000 total individuals, counting both targets and references;
a reference-assisted call may contain just one target. Each chromosome must
have at least two usable markers and a positive genetic span. Each sample
must have at least one observed call on each chromosome. Marker IDs must be
unique. The existing STvcf converter requires alphanumeric marker and
individual names. The map must be ordered by chromosome and nondecreasing
position. Common-stage markers must be monomorphic or biallelic; multiallelic
markers can be handled by the rare-marker stage, up to 100 declared alleles.
A selected multiallelic common marker is rejected. Small or sparse windows
can fail these conditions; use automatic windows or a larger explicit window.

\strong{Options and defaults.}
\describe{
  \item{\code{window = "auto"}}{Positive window length in cM, a vector with one
  entry per chromosome, \code{Inf} for whole-chromosome processing, or
  \code{"auto"} (the default) for the spacing rule below. Vector entries follow the
  chromosome order in the exported input; names are not used for matching.}
  \item{\code{overlap}}{Nonnegative overlap in cM; scalar or one value per
  chromosome. Automatic windows select overlap from marker spacing when it
  is omitted. Explicit windows use 2 cM when it is omitted. Finite windows
  must be at least 1.1 times their overlap.}
  \item{\code{window.markers = 4000000L}}{Integer cap on markers per finite
  window, at least 4. Explicit \code{window = Inf} ignores this cap. An
  automatic whole-chromosome choice exceeding the cap is rejected.}
  \item{\code{burnin = 3L}}{Number of burn-in iterations, an integer at least 1.}
  \item{\code{iterations = 12L}}{Number of subsequent phasing iterations,
  an integer at least 1.}
  \item{\code{phase.states = 280L}}{Maximum number of copying states;
  integer from 2 to 1000000. Increasing this can increase time and memory.}
  \item{\code{ne = 100000}}{Positive initial effective population size.
  EM may update the associated recombination parameter.}
  \item{\code{err = NULL}}{Initial allele mismatch probability. \code{NULL}
  selects the data-dependent initializer; an explicit value must lie in
  \eqn{[0,1]}. EM may update this value.}
  \item{\code{em = TRUE}}{Logical scalar controlling estimation of mismatch
  and recombination parameters.}
  \item{\code{initial.lr = 100000}}{Initial likelihood-ratio limit, at least 1.}
  \item{\code{step.scale = 3}}{Positive scale for donor-search step spacing.}
  \item{\code{rare = 0.002}}{Positive rare-allele threshold, at most 0.5.
  The algorithm uses sample carrier counts, with a minimum count boundary of
  three; this is not simply a filter on alternate-allele frequency.}
  \item{\code{buffer = 1}}{Positive donor-search context buffer in cM.}
  \item{\code{seed = NULL}}{Integer-valued seed in the range
  \eqn{[-2147483647,2147483647]}. \code{NULL} draws one seed from R's RNG,
  so \code{set.seed()} controls reproducibility. An explicit seed does not
  consume R random draws. The actual seed is returned in \code{settings}.}
}

\strong{Automatic windows.}
\code{window = "auto"} chooses context from the input map and the spacing of
common markers, without using phase truth. It targets a block of 400 common
markers. Fewer than 800 common markers on a chromosome leads to
whole-chromosome processing. Otherwise the window is at least 40 cM and spans
the longest local 400-marker interval. Automatic overlap is at least 2 cM and
spans the longest local 40-marker interval. A window that spans the chromosome
becomes \code{Inf}. Explicitly supplying \code{overlap} overrides the adaptive
overlap. These are context targets, not a guarantee of a minimum marker count
in every terminal window. Automatic windows are the default. They can use more memory than a
fixed 40-cM window when the marker map is sparse.

\strong{Reference requirements.}
Reference markers are matched by unique marker ID and reordered to target
order. Chromosome labels must agree and genetic positions must agree within
\eqn{10^{-6}} cM. The panel must contain complete calls at every target marker.
Homozygous panels have unambiguous phase. If any reference call is
heterozygous, the STvcf must have \code{attr(ref, "phased") = TRUE}; this is a
caller declaration, not independent verification of biological phase.

Declared allele labels are aligned, including reversed REF/ALT order. At least
one allele label must be shared at a matched marker. A reference-only allele
at an existing target marker can expand its allele definition. Reference-only
markers are ignored and reported. Reference individuals may have names also
used by targets in a separate data set. A loaded reference must be a different
data set from the target when writing back.

\emph{Allele coding matters:} SelectionTools' existing marker-data exporter
normalizes each marker's observed alleles to numeric labels and does not
retain original VCF nucleotide labels. Independently imported target and
reference panels must already use compatible exported allele coding. In
particular, independently normalizing a monomorphic panel can obscure which
allele it shares with a polymorphic target. Supplying an STvcf reference does
not recover discarded target labels. Check coding before phasing; this
interface does not infer strand orientation.

\strong{Missing data and write-back.}
The engine infers complete working haplotypes internally. With
\code{impute = FALSE}, originally missing calls are restored before return and
write-back. With \code{impute = TRUE}, inferred calls at these typed markers
are retained without a posterior-confidence cutoff.
No untyped markers, dosages, genotype probabilities or allele
probabilities are produced.

All computations precede write-back. The complete output is validated and
staged before the existing genotype arrays are changed. Original internal
allele codes, map and other dataset fields are retained. A failed phasing or
validation call does not partially write a chromosome into the source.
The returned STvcf uses exported allele coding. References are never written.
As with other genotype-changing operations, previously computed analyses
should be recomputed when their inputs have changed through imputation.

\strong{Threads and reproducibility.}
Use \code{st.set.num.threads} and \code{st.get.num.threads}.
There is no \code{nthreads} entry in \code{phasing.options}. OpenMP may grant
fewer workers than requested; the maximum actual team is reported. Builds
without OpenMP execute serially and report one phasing thread. OpenMP must
be enabled in the compiler and linker settings used to build the package. A fixed explicit seed gives deterministic
output across tested thread counts on the same build and input. Floating-point
and compiler differences can affect exact agreement across platforms.

}
\value{
Invisibly returns a list with:
\describe{
  \item{\code{data}}{Phased target-only STvcf, with \code{phased = TRUE}.
  Missing calls follow \code{impute}. This is also returned with
  \code{commit = FALSE}.}
  \item{\code{settings}}{Mode (\code{"cohort"} or \code{"reference"}),
  imputation policy, the supplied \code{ref}, effective
  \code{phasing.options}, requested
  \code{nthreads}, \code{actual.seed}, per-chromosome \code{window.settings}
  and \code{automatic.window}. Effective \code{ne} and \code{err} entries
  describe initialization controls, not final fitted EM estimates.}
  \item{\code{diagnostics}}{Total \code{windows}, \code{chromosome.windows},
  \code{maximum.worker.team}, \code{filled.genotypes} and
  \code{committed}. Reference runs additionally report reference sample and
  haplotype counts, invariant-check status, and ignored reference markers.}
}
}
\section{Subsetting}{
Restrict a copy of the source with \code{st.restrict.marker.data}
before calling \code{st.phase()} when only a subset of chromosomes or
individuals is wanted.
}

