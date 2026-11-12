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
  donor haplotypes. Reference haplotypes remain fixed. Further details appear in the reference
  requirements below. Default: \code{NULL}.}
  \item{phasing.options}{Uniquely named list of controls described below.
  Omitted entries use their defaults. Hyphens in option names are accepted as
  alternatives to dots, for example \code{"phase-states"}.
  Unknown or duplicate entries are errors.}
  \item{impute}{Logical scalar. \code{FALSE} preserves missing target
  genotypes. \code{TRUE} retains inferred hard calls at missing genotypes
  of existing target markers. Observed unordered genotypes are preserved
  under both policies. Default: \code{FALSE}.}
  \item{commit}{Logical scalar. \code{TRUE} writes the completed result to
  \code{data.set} and marks it phased. \code{FALSE} performs and validates
  phasing, discards the proposed genotypes, and returns only diagnostics without
  changing the data set. Default: \code{TRUE}.}
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
unique. The STvcf converter requires non-empty marker and individual names
without whitespace or control characters; ordinary punctuation is permitted. The map must be ordered by chromosome and nondecreasing
position. Common-stage markers must be monomorphic or biallelic; multiallelic
markers can be handled by the rare-marker stage, up to 100 declared alleles.
A selected multiallelic common marker is rejected. Small or sparse windows
can fail these conditions; automatic windows or a larger explicit window can satisfy the conditions.

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
  consume R random draws. The actual seed is reported with \code{commit = FALSE}.}
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
not recover discarded target labels. Coding compatibility is required before phasing; this
interface does not infer strand orientation.

\strong{Missing data and write-back.}
The engine infers complete working haplotypes internally. With
\code{impute = FALSE}, originally missing calls are retained in the validated proposal and
write-back. With \code{impute = TRUE}, inferred calls at these typed markers
are retained without a posterior-confidence cutoff.
No untyped markers, dosages, genotype probabilities or allele
probabilities are produced.

All computations precede write-back. The complete output is validated and
staged before the existing allele storage is changed. Original internal
allele codes, map, phenotypes and fitted effects are retained. Marker statistics
are refreshed with the committed genotypes; calculated design and relationship
matrices are invalidated. A failed phasing or
validation call does not partially write a chromosome into the source.
Loaded marker data are read directly in C into chromosome workspaces.
No complete STvcf snapshot, R genotype matrix or phased genotype object is
created as an internal interface or returned. References are never written.
As with other genotype-changing operations, previously computed analyses
should be recomputed when their inputs have changed through imputation.
Storage is selected automatically from the largest observed allele count at
any marker, using 4, 8 or 16 bits per diploid genotype plus dictionaries where
required. Imputation can promote storage when it introduces additional alleles.
The storage width is reported by \code{st.data.set.list()}.

\strong{Threads and reproducibility.}
\code{st.set.num.threads} and \code{st.get.num.threads} control and report the thread budget.
There is no \code{nthreads} entry in \code{phasing.options}. OpenMP may grant
fewer workers than requested; the maximum actual team is reported. Builds
without OpenMP execute serially and report one phasing thread. OpenMP must
be enabled in the compiler and linker settings used to build the package. A fixed explicit seed gives deterministic
output across tested thread counts on the same build and input. Floating-point
and compiler differences can affect exact agreement across platforms.

}
\value{
With \code{commit = TRUE}, invisibly returns \code{NULL} after updating the
supported data set. With \code{commit = FALSE}, returns a small diagnostics
list and discards the proposed genotypes. No genotype or full reference
object is returned. The dry-run diagnostics contain:
\describe{
  \item{\code{actual.seed}}{The actual seed used, including a seed drawn from R.}
  \item{\code{windows}, \code{chromosome.windows}}{Total and per-chromosome window counts.}
  \item{\code{maximum.worker.team}}{Maximum actual worker-team size.}
  \item{\code{requested.threads}}{Requested SelectionTools thread count.}
  \item{\code{filled.genotypes}}{Proposed number of filled genotypes.}
  \item{\code{committed}}{Always \code{FALSE} for returned dry-run diagnostics.}
  \item{\code{window.settings}}{Effective per-chromosome window choices.}
  \item{\code{reference}}{Reference sample and haplotype counts, invariant-check
  status and ignored-marker count, or \code{NULL} in cohort mode.}
}
}
\section{Subsetting}{
Restrict a copy of the source with \code{st.restrict.marker.data}
before calling \code{st.phase()} when only a subset of chromosomes or
individuals is wanted.
}

