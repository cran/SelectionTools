\name{st.phase.v1}
\alias{st.phase.v1}
\title{Statistical Phasing and Missing-Genotype Imputation of Diploid Marker Data}

\description{
Phases diploid biallelic SelectionTools marker data either from an external
SelectionTools reference data set or by reference-free cohort phasing.
Observed heterozygous genotypes receive a complete final phase;
posterior confidence is retained for diagnostics and for imputed genotypes,
but local ambiguity does not suppress otherwise valid chromosome phase.
Completely missing target genotypes can be imputed from posterior genotype
probabilities.  All proposed changes are constructed in private storage,
audited, and by default written back transactionally after a fresh final audit.
}

\usage{
st.phase.v1(
    reference.data.set = NULL,
    reference.individuals = NULL,
    initial.pair.weights = NULL,
    population.type = NULL,
    expected.breaks.per.morgan = NULL,
    max.states = NULL,
    error.rate = 0.001,
    break.rate1.per.morgan = NULL,
    break.rate2.per.morgan = NULL,
    confidence.threshold = 0.50,
    impute = TRUE,
    imputation.threshold = 0.90,
    individuals = NULL,
    chromosomes = NULL,
    burnin.max.iterations = 4,
    burnin.min.iterations = 2,
    burnin.convergence.fraction = 0.001,
    progressive.max.iterations = 2,
    progressive.min.iterations = 1,
    progressive.convergence.fraction = 0.001,
    burnin.states = 12,
    local.states = 8,
    stage1.minor.allele.count = 4,
    provisional.threshold = 0.80,
    commit = TRUE,
    data.set = "default"
)
}

\arguments{
  \item{reference.data.set}{
    \code{NULL} selects reference-free cohort phasing.  Otherwise, one
    non-empty character string naming an existing SelectionTools marker data
    set that supplies reference haplotypes.
  }

  \item{reference.individuals}{
    Individuals in \code{reference.data.set} to use as the reference panel.
    \code{NULL} selects all individuals.  Otherwise, supply unique individual
    names or unique one-based indices.  This argument is valid only in
    reference-data-set mode.

    All selected reference genotypes must be observed at every target marker.
    If every selected reference individual is homozygous, one haplotype is
    taken from each individual. If any selected reference individual is
    heterozygous, the two stored allele tracks of every selected reference are
    interpreted as already-phased haplotypes and \code{reference.data.set}
    must be marked phased. Data returned from the simulation backend and
    completely phased by \code{st.phase.v1()} carry this state automatically;
    consistently \code{|}-phased VCF input does as well.
  }

  \item{initial.pair.weights}{
    Optional nonnegative \eqn{R\times R} matrix of first-marker weights for
    ordered reference donor pairs.  It is valid only in reference mode and the
    sum of all entries must be positive.  When \code{NULL}, SelectionTools
    constructs a canonical upper-triangular labeling prior internally.  This
    breaks only the arbitrary global labels of the two copied homologues; it
    does not exclude an unordered donor pair.
  }

  \item{population.type}{
    Optional character string selecting a plant-population expectation for
    ancestry junctions.  Accepted values are \code{"DH"}, \code{"F2"},
    \code{"F3"}, \code{"F4"}, \code{"S2"}, \code{"S3"}, \code{"S4"},
    \code{"SSD"}, and \code{"Finfty"}, case-insensitively.

    The corresponding expected breaks per Morgan are 1 for DH and F2, 1.5 for
    F3 and S2, 2 for F4 and SSD, 1.75 for S3, and 1.875 for S4.
    \code{"Finfty"} denotes a general/outbred population and requests a
    lightweight data-derived copying-break estimate rather than a fixed plant
    pedigree value.

    If \code{population.type = NULL}, a reference panel consisting of exactly
    two selected homozygous individuals is recognized as \code{"F2"}; all
    other cases default to \code{"Finfty"}.
  }

  \item{expected.breaks.per.morgan}{
    Optional nonnegative numeric expert override for the expected ancestry or
    copying breaks per Morgan.  When supplied, this value overrides the value
    implied by \code{population.type} and is used for both internal break
    rates.  It cannot be combined with the legacy explicit pair
    \code{break.rate1.per.morgan}/\code{break.rate2.per.morgan}.
  }

  \item{max.states}{
    Positive whole number controlling the bounded donor set.  The default is
    \code{min(60,R)} in reference mode, where \eqn{R} is the number of
    reference haplotypes, and 60 in cohort mode.  In cohort mode it bounds the
    target-specific posterior donor subset and final imputation donor set.
  }

  \item{error.rate}{
    Allelic mismatch probability between zero and one.  The default is
    \code{0.001}.  For a biallelic marker, a copied donor allele emits the
    matching target allele with probability \eqn{1-e} and the alternative with
    probability \eqn{e}.  The parameter can absorb genotyping error, mutation,
    and imperfect representation by the donor panel.
  }

  \item{break.rate1.per.morgan, break.rate2.per.morgan}{
    Legacy advanced overrides for nonnegative copying-break intensities per
    Morgan.  They must be supplied together.  When either is supplied,
    \code{population.type} must not be supplied explicitly and
    \code{expected.breaks.per.morgan} must be \code{NULL}.

    For adjacent markers separated by \eqn{d} centiMorgans, intensity
    \eqn{\lambda} gives
    \deqn{q = 1 - \exp(-\lambda d/100)}
    as the probability of a copying break over the interval.  In reference
    mode the two values belong to the two ordered copying chains.  In cohort
    mode rate 1 is used for the Stage-1 scaffold and rate 2 for Stage-2
    low-frequency insertion.  When the legacy pair is not supplied, both are
    set to the effective value selected from \code{population.type},
    \code{expected.breaks.per.morgan}, or the \code{"Finfty"} estimator.
  }

  \item{confidence.threshold}{
    Probability threshold between 0.5 and 1 used for posterior-supported phase
    diagnostics and for assigning ordered phase to accepted imputed
    heterozygotes.  The default is \code{0.50}.  Numerically indistinguishable
    orientation probabilities are treated as a tie rather than allowing
    floating-point roundoff to create a hard phase call.  Observed
    heterozygotes are nevertheless completed from the continuous Viterbi path
    when an absolute local orientation is tied or falls below this threshold.
  }

  \item{impute}{
    Logical.  The default \code{TRUE} also decodes completely missing target
    genotypes from smoothed donor-state posteriors.  \code{FALSE} phases
    observed genotypes but leaves missing target genotypes missing.
  }

  \item{imputation.threshold}{
    Posterior genotype-probability threshold between zero and one for accepting
    a completely missing target genotype.  The default is \code{0.90}.
  }

  \item{individuals}{
    Target individuals to phase in reference mode.  \code{NULL} selects all;
    otherwise supply unique names or one-based indices.  Cohort mode is a
    whole-data-set procedure and requires \code{NULL}.
  }

  \item{chromosomes}{
    Target chromosomes to phase in reference mode.  \code{NULL} selects all;
    otherwise supply unique chromosome labels.  Cohort mode requires
    \code{NULL}.
  }

  \item{burnin.max.iterations}{
    Maximum number of synchronous Stage-1 Viterbi burn-in passes in cohort
    mode.  Default 4.  It must be at least \code{burnin.min.iterations}.
  }

  \item{burnin.min.iterations}{
    Minimum number of Stage-1 Viterbi burn-in passes before convergence can
    stop burn-in.  Default 2.
  }

  \item{burnin.convergence.fraction}{
    Cohort burn-in convergence tolerance.  After the minimum number of passes,
    burn-in stops when the fraction of changed consecutive-heterozygote phase
    relations is no greater than this value.  Default \code{0.001}.
  }

  \item{progressive.max.iterations}{
    Maximum number of synchronous Stage-1 posterior-refinement passes in
    cohort mode.  Default 2.
  }

  \item{progressive.min.iterations}{
    Minimum number of Stage-1 posterior-refinement passes before convergence
    can stop refinement.  Default 1.
  }

  \item{progressive.convergence.fraction}{
    Posterior-refinement convergence tolerance.  Default \code{0.001}.
  }

  \item{burnin.states}{
    Maximum number of balanced leave-one-out donor haplotypes retained per
    marker in the cohort Stage-1 Viterbi burn-in.  Default 12.
  }

  \item{local.states}{
    Number of local phase-invariant donor candidates retained per Stage-1
    marker while constructing the chromosome-wide target-specific donor
    ranking.  Default 8.
  }

  \item{stage1.minor.allele.count}{
    Minimum observed minor-allele count for a polymorphic biallelic marker to
    enter the iteratively phased Stage-1 scaffold.  Lower-frequency markers are
    inserted in Stage 2.  Default 4.
  }

  \item{provisional.threshold}{
    Probability threshold between 0.5 and 1 for allowing an allele at an
    originally missing genotype to become visible as a provisional donor allele
    in later cohort-refinement passes.  Default \code{0.80}.  This is distinct
    from final genotype imputation.
  }

  \item{commit}{
    Logical.  The default \code{TRUE} publishes the phased/imputed marker
    data transactionally after a fresh final audit.  Set \code{FALSE} for a
    staging-only run that performs the complete calculation and audit without
    modifying marker data.
  }

  \item{data.set}{
    One non-empty character string naming the SelectionTools marker data set to
    phase.  The default is \code{"default"}.
  }
}

\section{Default population models and break-rate selection}{
Normal end-user use does not require HMM break-rate tuning.  SelectionTools
first determines \code{population.type}.  Fixed plant-population types map to
expected ancestry-junction densities as follows:

\tabular{lrl}{
\strong{Type} \tab \strong{Breaks/Morgan} \tab \strong{Interpretation}\cr
DH  \tab 1.000 \tab doubled haploids directly from an F1\cr
F2  \tab 1.000 \tab F2 from an F1\cr
F3  \tab 1.500 \tab one additional random-mating generation\cr
F4  \tab 2.000 \tab two additional random-mating generations\cr
S2  \tab 1.500 \tab two generations of selfing from the F1\cr
S3  \tab 1.750 \tab three generations of selfing\cr
S4  \tab 1.875 \tab four generations of selfing\cr
SSD \tab 2.000 \tab single-seed descent / near-inbred selfing limit\cr
Finfty \tab estimated \tab general or outbred population\cr
}

For \code{"Finfty"}, SelectionTools estimates a single copying-break intensity
on the same per-Morgan scale as the Li--Stephens transition.  The estimator is
lightweight and deterministic: it uses at most 500 individuals, common
biallelic markers, and a bounded sample of short- to intermediate-distance
marker pairs, then fits a one-dimensional LD-decay curve.  It does not rerun
phasing and does not perform a multi-parameter grid search.  The resulting
value is returned as \code{expected.breaks.per.morgan} and copied to both
internal break rates unless the legacy pair is explicitly supplied.

Precedence is: the explicit legacy two-rate pair, otherwise
\code{expected.breaks.per.morgan}, otherwise the fixed value implied by
\code{population.type}, with \code{"Finfty"} triggering estimation.
}

\section{Phasing modes}{
A named \code{reference.data.set} invokes the ordered donor-pair reference HMM.
Reference and target markers are matched by marker name.  With
\code{reference.data.set = NULL}, the target cohort supplies leave-one-out
provisional donor haplotypes and is phased with the two-stage common-marker / low-
frequency procedure.

The cohort-only iteration/state arguments are ignored by the C reference
pipeline and must not be supplied explicitly in reference mode.  In cohort
mode \code{reference.individuals}, \code{initial.pair.weights},
\code{individuals}, and \code{chromosomes} are not accepted.
}

\section{Scope and restrictions}{
The public phaser operates on diploid SelectionTools marker data with at most
two biological alleles per marker.  Completely missing genotypes must be
missing on both stored alleles; partial diploid missingness is rejected.  Map
positions are in centiMorgans and chromosomes are processed separately in map
order.

Reference genotypes must be complete at all target markers.  Heterozygous
reference individuals are accepted only as already-phased diploids.  A
reference allele absent from the corresponding target marker's biological
allele dictionary is rejected rather than silently recoded.

Cohort initialization can use private provisional alleles at completely missing
genotypes.  Such provisional alleles are never interpreted as observed
biological genotypes and are never directly committed as imputation calls.
}

\section{Statistical theory and implemented algorithm}{

\subsection{Notation and copying model}{
For chromosome markers \eqn{m=1,\ldots,M} at positions \eqn{x_m} centiMorgans,
let \eqn{g_m=\{a_m,b_m\}} be the unordered target genotype and \eqn{H_{mr}}
the allele carried by donor haplotype \eqn{r}.  Biological SelectionTools
alleles are translated to compact per-marker states for phasing and translated
back only during staging.

For allelic mismatch probability \eqn{e}, a copied donor allele emits its match
with probability \eqn{1-e}; at a biallelic mismatch the probability is
\eqn{e}.  A completely missing observation is neutral.  For interval length
\eqn{\Delta_m=(x_m-x_{m-1})/100} Morgans and copying-break intensity
\eqn{\lambda},
\deqn{q_m=1-\exp(-\lambda\Delta_m).}
A retained donor persists with mass \eqn{1-q_m} and the break mass \eqn{q_m}
is redistributed over the current bounded donor states.
}

\subsection{Reference-panel donor selection and ordered diploid HMM}{
Reference donor selection is invariant to the arbitrary input order of a target
heterozygote.  Contiguous compatibility runs between each donor haplotype and
either allele of the unordered target genotype provide an IBS ranking; at most
\code{max.states} distinct donors are retained per marker.

The hidden diploid state is the ordered donor pair \eqn{Z_m=(r_m,s_m)}.  The
two donor chains transition independently.  For a heterozygote,
\deqn{
E_m(r,s)=P(a_m\mid H_{mr})P(b_m\mid H_{ms})+
         P(b_m\mid H_{mr})P(a_m\mid H_{ms}).
}
Forward/backward posterior probabilities and a Viterbi donor-pair path are
computed with factorized \eqn{O(MK^2)} transitions after candidate construction.
}

\subsection{complete reference phase}{
Posterior donor-pair mass is first partitioned between the two possible ordered
orientations of each observed heterozygote.  Calls that reach
\code{confidence.threshold} are retained only when the two posterior
orientation probabilities are numerically distinguishable.  This prevents
roundoff around an exact label-symmetric 50:50 posterior from creating an
artificial hard orientation.

Observed heterozygotes that remain locally unresolved are then completed rather
than omitted.  The deterministic Viterbi donor-pair traceback provides a
continuous labeled chromosome path.  At an unresolved marker the two ordered
emission components for the Viterbi state are compared; if one is larger it
sets the orientation, and an exact tie retains the running chromosome
orientation.  This continue/swap rule avoids
mistaking arbitrary homologue-label symmetry for absence of relative phase.
Thus ordinary observed heterozygotes receive complete phase output.
}

\subsection{Phased heterozygous reference individuals}{
When all selected reference individuals are homozygous, one haplotype per
individual is used. If at least one selected reference is heterozygous, the
reference marker data set must have its phased state set to \code{TRUE}
and both stored allele tracks of every selected reference individual are
used as phased homologues. Hence \eqn{N_R} selected diploids contribute
\eqn{2N_R} donor haplotypes in that case.
}

\subsection{PBWT initialization of a reference-free cohort}{
Reference-free mode starts from unordered genotypes and never treats the input
allele-1/allele-2 order as known phase.  A PBWT-based phase-invariant
initialization scores candidate allele assignments from compatible prefix
matches to other individuals.  Completely missing genotypes receive private
provisional major-allele proxies so initialization can continue; the original
missingness mask is retained.
}

\subsection{Two-stage cohort phasing}{
Polymorphic markers with observed minor-allele count at least
\code{stage1.minor.allele.count} form Stage 1.  Rarer polymorphic markers are
withheld and inserted once in Stage 2.  Every iterative cohort pass is
leave-one-individual-out and synchronous: a pass reads a frozen snapshot and
writes a separate next snapshot, so individual processing order cannot alter
that pass.
}

\subsection{Stage-1 burn-in and posterior refinement}{
Stage 1 begins with synchronous Viterbi burn-in using at most
\code{burnin.states} balanced phase-invariant donor states.  Burn-in stops after
at least \code{burnin.min.iterations} when its changed-link fraction is no
greater than \code{burnin.convergence.fraction}, or at
\code{burnin.max.iterations}.

Posterior refinement then retains \code{local.states} phase-invariant local
candidates per marker, aggregates support across the chromosome, and keeps at
most \code{max.states} target-specific donors.  Originally missing donor alleles
can enter later frozen panels only when their haploid posterior reaches
\code{provisional.threshold}.  Refinement uses analogous minimum, maximum, and
convergence controls.
}

\subsection{Stage-2 low-frequency insertion and ties}{
Stage-2 markers are oriented once conditional on the final common-marker
scaffold.  Haploid donor-state probabilities are computed at neighboring
Stage-1 markers and interpolated along genetic distance.  For a heterozygous
low-frequency genotype, the two ordered allele products are compared.

If one orientation has larger support it is chosen.  If the two scores are
exactly equal, including a singleton minor allele absent from every
leave-one-out donor, a stable 50:50 pseudo-random tie decision based on the
individual, marker, and chromosome coordinates chooses one orientation.  This
resolves Stage-2 ties while preserving
bitwise reproducibility across OpenMP thread counts.
}

\subsection{Complete cohort output}{
The Stage-1 relative-phase chain and Stage-2 insertion decisions are always
written to the private stage for observed heterozygotes.  Confidence flags are
retained as diagnostics but do not veto an otherwise valid chromosome or turn
a local singleton ambiguity into chromosome-wide missing phase.  Consequently
normal cohort phasing has complete observed-heterozygote phase coverage.
}

\subsection{Posterior missing-genotype imputation}{
When \code{impute = TRUE}, final biological imputation is separate from the
provisional missing-data bootstrap.  In cohort mode the target is excluded and
only donor individuals whose biological genotype was originally observed at
every marker on that chromosome are eligible for the imputation HMM.  If no
such donor remains, missing genotypes stay unresolved rather than causing the
pipeline to fail.

The smoothed ordered-state posterior is collapsed to probabilities for
unordered diploid genotypes.  A genotype is accepted only when its largest
posterior reaches \code{imputation.threshold}.  For an accepted imputed
heterozygote, ordered phase is attached only when the conditional orientation
posterior reaches \code{confidence.threshold}; otherwise the genotype can be
imputed without claiming confident ordered phase.
}

\subsection{Transactional staging, audit, and commit}{
Phasing and imputation operate on private storage.  The audit verifies data
sizes, marker mapping, map positions, source-byte freshness, preservation of
observed unordered biological genotypes, and consistency of staged writes.
With \code{commit = TRUE}, a fresh audit is performed immediately before
publication. With \code{commit = FALSE}, the complete calculation and audit
run without changing the source marker data. A successful complete-data-set
commit with no unresolved heterozygous phase marks the SelectionTools marker
data set as phased. Partial phasing requests do not by themselves promote an
otherwise unphased data set to the phased state.
}

\subsection{Computational structure and OpenMP}{
Reference diploid HMM work is \eqn{O(MK^2)} after bounded donor selection.
Cohort mode repeats target-specific leave-one-out haploid calculations across a
small fixed number of synchronous burn-in/refinement passes and inserts rare
markers once.  When OpenMP is available, independent target individuals are
parallelized within frozen passes and Stage-2 insertion.  Chromosomes are not
parallelized.
}
}

\section{Interpretation of defaults and expert controls}{
The default interface is intended for end users rather than per-data-set
parameter tuning.  Known plant population types use fixed biologically
interpretable ancestry-junction expectations.  General populations use the
lightweight \code{"Finfty"} estimator.  The default \code{error.rate = 0.001},
state bounds, and short convergence-controlled iteration limits are fixed.

\code{expected.breaks.per.morgan} is the preferred expert override when the
expected ancestry-junction density is known.  The two legacy break-rate
arguments remain available only for advanced/backward-compatible analyses that
need distinct chain or stage rates.  Increasing state bounds generally uses
more donor information at higher computational cost.  Increasing
\code{imputation.threshold} reduces the number of missing genotypes accepted.
}

\value{
An object of class \code{"st.phase.result"} and \code{"list"}.  Full
marker-wise posterior arrays are not returned.  Common diagnostics include:
\describe{
  \item{\code{success}}{Logical indicating successful completion.}
  \item{\code{mode}}{\code{"reference"} or \code{"cohort"}.}
  \item{\code{pipeline.status}, \code{pipeline.status.name}}{Pipeline status.}
  \item{\code{audit.status}, \code{audit.status.name}}{Private-stage audit status.}
  \item{\code{commit.requested}, \code{commit.attempted},
    \code{commit.status}, \code{commit.status.name}, \code{committed}}{Commit diagnostics.}
  \item{\code{n.ind}, \code{n.marker}, \code{n.chrom}, \code{n.genotype}}{Target dimensions.}
  \item{\code{n.missing}, \code{n.homozygous}, \code{n.heterozygous}}{Source genotype counts.}
  \item{\code{n.written}, \code{n.changed.genotypes}}{Staged phase/imputation write counts.}
  \item{\code{n.imputation.missing.requested}, \code{n.imputed},
    \code{n.imputation.unresolved}}{Missing-genotype imputation counts.}
  \item{\code{population.type}}{Effective population type.}
  \item{\code{expected.breaks.per.morgan}}{Effective single break expectation;
    \code{NA} when the legacy two-rate pair was supplied.}
  \item{\code{breaks.estimated}}{Logical indicating whether \code{"Finfty"}
    estimation supplied the effective single break rate.}
  \item{\code{break.rate1.per.morgan}, \code{break.rate2.per.morgan}}{Effective
    break rates passed to the C phaser.}
  \item{\code{error.rate}, \code{max.states}}{Effective mismatch rate and donor-state bound.}
}

Reference mode additionally reports reference-panel sizes and reference-mode
run diagnostics.  Cohort mode additionally reports Stage-1/Stage-2 marker
counts, burn-in/refinement iteration diagnostics, provisional-state counts,
relative same/switch-link diagnostics, and confidence diagnostics.  Under the
complete-output policy, low posterior confidence for an observed
heterozygote is diagnostic rather than an instruction to leave that observed
phase unwritten.
}

