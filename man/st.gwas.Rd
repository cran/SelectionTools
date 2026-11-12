\name{st.gwas}
\alias{st.gwas}
\alias{print.st.gwas}
\alias{as.data.frame.st.gwas}
\alias{plot.st.gwas}
\alias{compare.gwas.fits}
\alias{summarize.gwas.fit}
\alias{summarize.gwas.simulation}
\alias{print.st.gwas.summary}
\title{Single-Marker Association Scans with Mixed Models}
\description{Tests additive marker effects with optional genomic
relationships, population covariates and known residual covariance.
Optionally calibrates a genome-wide significance threshold by parametric
simulation from the fitted null model.}
\usage{
st.gwas(fixed, random, data, marker.data.set = "default",
        K = NULL, n.PC = 0L, min.MAF = 0.05, P3D = TRUE,
        max.missing = 1, missing = "omit", map = NULL, rcov = NULL,
        control = list(), sim.threshold = list())
compare.gwas.fits(..., adjust = "FDR", alpha = 0.05, trait = NULL)
summarize.gwas.fit(fit, adjust = "FDR", alpha = 0.05,
                   top = 6L, trait = NULL)
summarize.gwas.simulation(fit, alpha = 0.05, top = 6L, trait = NULL)
\method{print}{st.gwas.summary}(x, digits = NULL, char.per.line = NULL, ...)
\method{as.data.frame}{st.gwas}(x, ...)
\method{print}{st.gwas}(x, digits = NULL, top = NULL, bottom = NULL,
                      char.per.line = NULL, ...)
\method{plot}{st.gwas}(x, which = "manhattan", trait = NULL,
                     alpha = 0.05, adjust = "BH", threshold = NULL,
                     chromosome = NULL, xlim = NULL,
                     position.scale = 1e6, position.unit = "Mb",
                     shade = NULL, highlight = NULL, labels = FALSE,
                     level = NULL, main = NULL,
                     col = c("grey35", "cornflowerblue"),
                     cex = 0.5, pch = 16, ...)
}
\arguments{
\item{fixed}{Two-sided fixed-effects formula, including numeric covariates
and factors as needed. A \code{cbind} response requests separate univariate
scans for its columns, not a joint multivariate test.}
\item{random}{Required one-sided formula identifying the genotype column,
for example \code{random = ~gen}. Exactly one untransformed column name
is supported. The polygenic effects for these genotype levels have
covariance proportional to \code{K}; \code{K = "A"} constructs that
relationship internally from the marker panel. The column contains
nonmissing marker-individual identifiers. Repeated records for a genotype
share its polygenic effect. With \code{K = FALSE}, the polygenic effect
is omitted, but the formula still identifies marker individuals.}
\item{data}{Phenotype data frame. Repeated records per individual are
allowed. Column names beginning with \code{.st.gwas.} are reserved.}
\item{marker.data.set}{compiled marker dataset name, or numeric matrix with
individuals in rows and markers in columns. Matrix row and column names
must be unique identifiers. Values are allele dosages between zero and two,
including fractional imputed dosages; \code{NA} denotes missing calls.
Converts rrBLUP scores in minus one, zero, one coding by adding one.
compiled datasets are read directly without exporting the complete matrix.}
\item{K}{NULL (the default), or FALSE, omits kinship adjustment. The required
\code{random} formula still identifies marker individuals; no polygenic
random effect is fitted when kinship is omitted. The quoted string \code{"A"}
builds a VanRaden relationship internally from \code{marker.data.set}, using
all polymorphic biallelic markers with mean imputation. Symbolic \code{K = "A"}
does not evaluate a caller object named \code{A}. Unquoted \code{K=A}
evaluates the actual R object \code{A}. An alternative is a
named finite symmetric positive-semidefinite matrix covering every marker
individual. Its axes are aligned by identifier. Kinship accumulation,
eigendecomposition and PC calculation use compiled buffers. Exported copies
are retained in the returned null models for compatibility.}
\item{n.PC}{Number of leading eigenvectors of \code{K} added as fixed
population covariates. Requires explicit \code{K = "A"} or a supplied relationship matrix. For an ordinary
scan with supplied population covariates, \code{fixed} contains those covariates
and \code{K=FALSE,n.PC=0} specifies the scan.}
\item{min.MAF}{Minimum minor allele frequency, between zero and one half.
Frequency is computed over unique phenotyped individuals retained for the
trait and available fixed predictors, not over replicated records.}
\item{P3D}{TRUE fits the null covariance once and holds its shape fixed
during the scan. FALSE re-estimates free covariance parameters per marker.
The same setting is used for every simulated scan. Each simulated response
receives its own null-model fit with \code{method = "parametric"}.
\code{method = "parametricFixedVC"} requires P3D = TRUE and reuses the
original null covariance across replicates.}
\item{max.missing}{Maximum missing genotype fraction per marker, over
unique retained phenotyped individuals.}
\item{missing}{\code{"omit"} removes all records of individuals with a
missing call at the tested marker. \code{"mean"} explicitly imputes their
dosage by the observed mean for retained individuals.}
\item{map}{NULL, or a data frame with \code{marker}, \code{chromosome} and
\code{position} columns matched by marker identifier. A compiled dataset
supplies its stored physical map by default, falling back to its genetic
map when no physical map is stored. Physical coordinates imported from
the STvcf \code{physical.map} attribute need no separate \code{map} argument. Unmapped markers use scan order in
plots. Positions retain the units of the supplied map.}
\item{rcov}{Optional known residual covariance matrix, on its absolute
response scale. Naming, alignment and restrictions are those of
\code{\link{st.lmm}}. No additional residual variance is fitted.}
\item{control}{Fitting controls passed to the compiled model. Additionally
\code{batch.size} sets the number of marker columns processed together
(default 64; from 1 to 4096).}
\item{sim.threshold}{A list controlling simulated significance thresholds.
The default empty list disables simulation and does not use the random
number generator. A nonempty list enables simulation and accepts only
the uniquely named entries \code{nsim}, \code{alpha}, \code{seed}, and
\code{method}. The method is either \code{"parametric"} (the default) or
\code{"parametricFixedVC"}.
Omitted entries default to \code{nsim = 999L}, \code{alpha = 0.05}, and
\code{seed = NULL}. \code{nsim} must be a positive integer and \code{alpha}
must be a finite scalar strictly between zero and one. The number of
replicates must satisfy \code{floor(alpha * (nsim + 1)) >= 1}.
\code{nsim} must be smaller than \code{.Machine$integer.max}, and a
non-NULL \code{seed} must be a finite integer with absolute value no
larger than \code{.Machine$integer.max}. An explicit seed restores the
caller's previous random-number-generator kind and state, including the
absence of \code{.Random.seed}. With NULL, one master seed is drawn from
the current generator and recorded; only that draw advances the caller's
state. Scoped simulation does not support Box--Muller or user-supplied
random-number generators. Unknown, duplicated, or unnamed list entries
are errors.}
\item{fit}{One fitted \code{st.gwas} object to summarize.}
\item{x}{An \code{st.gwas} object, or a \code{st.gwas.summary} result for its print method.}
\item{digits,top,bottom,char.per.line}{For summaries, \code{top} is a positive
integer number of leading markers, supplemented by any significant chromosome
leads outside this set. Otherwise these are compiled table printing controls;
row limits may be \code{Inf}.}
\item{which}{One or more of \code{"manhattan"}, \code{"qq"}, and
\code{"region"}, in the requested order. Plotting does not change the
layout of the graphics device. A regional plot requires \code{chromosome}.}
\item{trait}{NULL for all traits, or one trait name to plot or summarize.}
\item{alpha}{Significance level for summary decisions or a BH, Bonferroni, or unadjusted plotted
threshold. When plotting a simulated threshold, an omitted \code{alpha}
uses its stored level. An explicitly supplied different level is an error;
recalibration is requested through \code{sim.threshold} in \code{st.gwas}.}
\item{adjust}{For summaries, \code{"FDR"} (the default), its alias \code{"BH"},
\code{"none"}, or \code{"simulation"}. FDR uses the stored trait-wise
Benjamini--Hochberg probabilities. Summary output has no Bonferroni column.
For plotting: \code{"BH"}, \code{"bonferroni"}, \code{"none"}, or
\code{"simulation"} for Manhattan and regional thresholds. When this
argument is omitted and simulation results are present, the stored
simulated threshold is selected automatically. Without simulation results,
the default is BH. Supplying an adjustment explicitly overrides this
automatic selection. BH is the Benjamini--Hochberg procedure; rrBLUP plots
use Storey q-values. An explicit \code{"simulation"} request requires a
complete simulated threshold for each requested trait.}
\item{threshold}{NULL uses the automatic or explicitly selected adjustment;
FALSE omits the horizontal significance line. Alternatively, one
or more raw p-values in (0, 1] specify the drawn thresholds. FALSE and numeric
values override \code{adjust}, including automatic simulated thresholds.
The stored simulated cutoff is drawn directly on its minus-log10 scale,
so its raw p-value may underflow to zero without losing the line.
Plotting never runs or repeats simulations.}
\item{chromosome}{One chromosome identifier for a regional plot.}
\item{xlim}{Two limits for a regional plot, in the original units of the
supplied map. The threshold is still computed from the complete trait scan.}
\item{position.scale}{Positive divisor converting map positions to regional
plot coordinates. The default \code{1e6} converts base pairs to megabases.}
\item{position.unit}{Unit label for the regional horizontal axis. This label
and \code{position.scale} must agree with the supplied map.}
\item{shade}{NULL, or two positions in the original map units delimiting a
shaded interval in a regional plot.}
\item{highlight}{NULL, or marker identifiers to highlight in coral.}
\item{labels}{TRUE labels the highlighted markers; FALSE omits labels.
A character vector labels the specified marker identifiers.}
\item{level}{NULL omits a QQ envelope. A probability between zero and one
adds a pointwise beta-order-statistic reference envelope. This envelope
assumes independent uniform null p-values; marker linkage disequilibrium
means that it is not a simultaneous genome-wide confidence band.}
\item{main}{NULL for the trait name, or a plot title.}
\item{col}{Plot colours. A Manhattan plot cycles these colours by
chromosome; QQ and regional plots use the first colour.}
\item{cex,pch}{Point size and plotting symbol for customized plots.}
\item{...}{For \code{compare.gwas.fits}, individual fitted \code{st.gwas}
objects. Named arguments supply model labels; otherwise bare fitted-object
names are used, with fit1, fit2, etc. as fallback labels for expressions.
Model labels must be unique; name repeated inputs explicitly.
For customized plots, named graphical arguments passed to
\code{graphics::plot}, such as \code{xlab}, \code{ylab}, \code{ylim},
and \code{bty}; \code{x}, \code{y}, and \code{type} are reserved.
Reserved for the printing and data-frame methods.}
}
\details{
Without \code{rcov}, tests use an F statistic with numerator df one and
denominator df equal to the marker model residual df. With P3D the null
covariance shape is retained and the marker model residual scale is
re-estimated. This follows the testing convention of rrBLUP's GWAS
function for complete biallelic genotypes. Known absolute \code{rcov}
instead gives a one-df GLS Wald chi-squared test, without rescaling the
entire covariance matrix.

For omitted genotypes, the covariance is restricted to retained
observations and solved again. A principal submatrix of the original
inverse is not substituted for this calculation. Thus missing-call
results can differ from rrBLUP's P3D shortcut. Fixed-effect rank and
estimability are checked for each tested marker.

Multiallelic, monomorphic, missing, filtered or nonestimable loci have a
status and no association statistic. Marker-specific fitting failures
produce a message and allow the remaining scan to continue. compiled
reference/effect allele values are stored biological allele codes, not
necessarily nucleotide strings. With strictACTG input they use
A=1, C=2, G=3, T=4. Numeric-matrix effects are per unit of supplied dosage.

BH and Bonferroni adjustments are computed separately per trait over
successful tests. Bounded batches avoid an additional complete dosage
export, but fitting still requires the dense relationship/covariance
storage used by the selected solver. The global thread setting from
\code{\link{st.set.num.threads}} is respected and is not changed.
The scan does not change marker data or stored prediction state.
\code{st.gwas} returns its result without opening a plot; plotting is a
separate call to \code{plot}.

\subsection{Simulated significance thresholds}{
With a nonempty \code{sim.threshold} list, the original marker-free null
model generates independent Gaussian response replicates. Its fitted
fixed mean and covariance parameters define the generating model. Random
effects and errors are regenerated jointly, using the unconditional
simulation convention \code{re.form = NA} of \code{\link{simulate.st.lmm}}.
The realized fitted random effects are not held fixed.

The marker data, resolved relationship matrix, principal components,
fixed-effects design, phenotype-to-individual mapping, and marker filters
remain those of the original analysis. Original missing response cells
are restored before each scan. Responses expressed through a formula
transformation are simulated on the retained analysis scale, with no
second application of that transformation. Repeated observations retain
their common individual mapping and its covariance structure. Known
absolute \code{rcov} retains the GLS Wald testing convention described
above; otherwise the simulated scans use the same F tests as the original.

With \code{method = "parametric"}, each replicate refits its null model.
With \code{P3D = TRUE}, its fitted
covariance shape is reused over that replicate's markers and each marker
retains its residual-scale estimate. With \code{P3D = FALSE}, free
covariance parameters are refitted per marker. Missing genotype calls
are handled by the original \code{missing} setting, including the
restricted-covariance solve for omitted observations.

With \code{method = "parametricFixedVC"}, the original fitted variance
components are held fixed for every replicate and marker; no simulated
model undergoes variance-component estimation. This method requires
\code{P3D = TRUE}. The null covariance factorization and complete or
mean-imputed marker information are reused. Fixed effects are projected
out separately for each response, and each marker retains its
response-specific residual-scale estimate for F tests. Omitted genotype
calls use a restricted-covariance solve with the same original variance
components. This is a faster conditional, plug-in calibration; it does
not reproduce the variation from re-estimating variance components.
Both methods generate Gaussian responses from the same original fitted
null model and use the same marker family and Monte Carlo calculation.

For each trait, the eligible marker family is fixed by the original
analysis and the maximum compiled minus-log10 p-value is recorded for every
simulated genome scan. Here \eqn{B} is \code{nsim}, \eqn{M_b} denotes these
maxima, and \eqn{s} is an observed marker score. The simulation-adjusted p-value is
\deqn{p_{\mathrm{sim}}(s) =
       \frac{1 + \sum_{b=1}^{B} I(M_b \geq s)}{B + 1}.}
The comparison includes ties. The definitions are
\eqn{h = \lfloor \alpha(B+1) \rfloor} and \eqn{k = B+1-h}.
The score cutoff is the \eqn{k}th ordered maximum in ascending order;
scores strictly above it meet the requested level. Significance decisions
use \code{p.sim <= alpha}. The cutoff is an observed order statistic,
without interpolating between maxima. Its raw p-value is supplied for
convenience, but the score remains authoritative if the raw value
underflows.

The nominal level concerns the maximum over markers within one trait
under its fitted global null model. Fitted generating parameters and a
finite simulation count make this a model-based approximation. Every
trait is calibrated separately; the level does not jointly cover several
traits or different model choices.

Failures of eligible marker tests or replicate fits, including nonfinite
simulated scores, are retained in the diagnostics. Failed draws are not
replaced, and failed markers are not
silently removed from the simulated maximum. If a trait's calibration is
incomplete, its threshold and simulation-adjusted p-values are unavailable.
The maxima, settings, and failure diagnostics remain in the returned
object. A trait with no eligible markers has no simulated threshold.
Only replicate maxima are retained, not every simulated marker table.
}

The existing Manhattan and QQ plotting calls are retained. Customized
Manhattan, QQ and regional plots use standard R graphics without additional
packages. A mapped Manhattan plot orders numeric chromosome labels
numerically and keeps the encounter order of other chromosome labels.
If any successful marker lacks a map position or chromosome, it uses
marker order instead. Regional plots require mapped markers. The BH line
identifies the largest raw p-value among BH discoveries, if any; it
represents false-discovery-rate control and is not a familywise threshold.
The Bonferroni line uses the number of successful tests in the complete
trait scan, even when plotting a region.

A completed simulated threshold is used automatically for Manhattan and
regional plots unless \code{adjust} or \code{threshold} supplies an
override. It is shown as a coral dashed line labelled with the simulation
level and calculated score cutoff. The automatic vertical range includes
that line. Regional significance colouring uses the simulation-adjusted
p-values and the stored level. If automatic selection encounters an
incomplete calibration, a warning is given and the significance line is
omitted. An explicit \code{adjust = "simulation"} request instead fails
with an error. Neither route substitutes a BH threshold for an incomplete
simulated threshold.

For customized plots, exact zero raw p-values are clipped to the smallest
positive normal double when drawing QQ points, explicitly supplied
raw-p-value thresholds, or association points with nonfinite scores.
Finite compiled association scores and stored simulated score cutoffs are
used directly.
}
\section{compiled GWAS summaries}{
\code{compare.gwas.fits} summarizes successful finite marker tests separately
for each model and trait. It returns the lead (smallest raw p-value) on each
chromosome containing markers significant at \code{alpha} under the requested
adjustment. This is a chromosome-level report, not a procedure for separating
multiple linked peaks on one chromosome. Without discoveries it reports the
overall lead; a trait with no successful tests has counts but no marker rows.
\code{summarize.gwas.fit} returns the first \code{top} successful markers
ordered by raw p-value, retaining additional significant chromosome leads.
Ties retain original result-row order. No fit is refitted.
\code{summarize.gwas.simulation} uses \code{p.sim} and rejects absent or
incomplete calibration for any selected trait. It requires status ok, all
requested simulations successful, zero failures and valid adjusted probabilities.
All counting, ordering, selection and validation use compiled C workspace.
}
\value{The three summary functions return a serializable
\code{st.gwas.summary} list with \code{models} (model, trait, PCs, tested,
significant), \code{markers} (model, trait, chromosome, marker, position,
MAF, effect, p, p.adjusted, score), \code{adjust}, \code{alpha} and
\code{comparison}. Printed adjusted probabilities are labelled FDR or p.sim
as appropriate. Printing uses the common compiled table printer, respects
\code{char.per.line} and prints all returned rows. Probabilities retain
\code{digits} decimal places, with values below the display limit shown
as, for example, \code{<0.0001} at four decimal places.

A serializable \code{st.gwas} list containing \code{results},
\code{null.models}, trait labels and scan settings. The results table has
\code{trait}, \code{marker}, \code{chromosome}, \code{position},
\code{reference.allele}, \code{effect.allele}, \code{MAF}, \code{n.lines},
\code{nobs}, \code{n.missing}, \code{effect}, \code{SE}, \code{statistic},
\code{DF}, \code{p}, \code{score}, \code{p.BH},
\code{p.bonferroni}, \code{status} and \code{message}. Here \code{score}
is minus log10 of the raw p-value and \code{DF} is the denominator DF
(infinite for the Wald test).

When simulation is enabled, \code{results} additionally contains
\code{p.sim}, and \code{simulation} contains:
\describe{
\item{\code{method, nsim, alpha, seed}}{The simulation method, requested
replicate count, nominal level, and recorded master seed.}
\item{\code{thresholds}}{A data frame with columns \code{trait},
\code{alpha}, \code{nsim}, \code{successful}, \code{failed},
\code{n.markers}, \code{score}, \code{p}, and \code{status}.
Status is \code{"ok"}, \code{"incomplete"}, or \code{"no.markers"}.
\code{score} is the simulated cutoff and \code{p} its raw-p-value
equivalent.}
\item{\code{max.score}}{A numeric matrix with \code{nsim} rows and one
named column per trait. Failed replicates have missing maxima.}
\item{\code{diagnostics}}{A data frame identifying the trait, one-based
replicate number, status, message, and marker if a marker test failed.}
}
Without simulation, the ordinary result layout is retained and no
\code{simulation} component is added.
\code{as.data.frame} returns the marker results table; printing and
plotting return the object invisibly.}
\references{
Endelman JB (2011). Ridge regression and other kernels for genomic selection
with R package rrBLUP. The Plant Genome 4:250--255.

Yu J et al. (2006). A unified mixed-model method for association mapping
that accounts for multiple levels of relatedness. Nature Genetics 38:203--208.

Kang HM et al. (2010). Variance component model to account for sample
structure in genome-wide association studies. Nature Genetics 42:348--354.

Joo JWJ, Hormozdiari F, Han B, Eskin E (2016). Multiple testing correction
in linear mixed models. Genome Biology 17:62.
\doi{10.1186/s13059-016-0903-6}

Phipson B, Smyth GK (2010). Permutation P-values should never be zero:
calculating exact P-values when permutations are randomly drawn.
Statistical Applications in Genetics and Molecular Biology 9:Article 39.
\doi{10.2202/1544-6115.1585}
}

