\name{propsd}
\alias{propsd}
\alias{print.propsd}
\title{Proportional Loss of Expected Sample Standard Deviation}
\description{
Evaluates the proportional loss of expected sample standard deviation under
a fixed Gaussian covariance model. Exact sample-variance moments determine
a two-moment Gamma approximation to the expected sample standard deviation.
}
\usage{
propsd(G, parents = NULL, counts = NULL)

\method{print}{propsd}(x, digits = NULL, ...)
}
\arguments{
\item{G}{Complete finite symmetric positive-semidefinite covariance-proportional
matrix with at least two reference parents and positive centered contrast
variance. Singular matrices are supported. Both dimensions use the same order.
Row and column names, when present, are complete, unique and match in order.
The reference contains one occurrence per row. Its covariance scale and
allele-frequency reference remain fixed.}
\item{parents}{Selected row indices or parental names, with repetitions
permitted. Repeated occurrences represent the same parental random variable.
At least two selected occurrences are required. Cannot be combined with
\code{counts}. Character labels for an unnamed matrix are \code{"1"},
\code{"2"}, and so on.}
\item{counts}{Alternative nonnegative integer multiplicities in full
reference order. Optional names match that order exactly. Each count and their
sum must be at most \code{2^53 - 1}; the sum is at least two.
Fractional gene contributions are not multiplicities.}
\item{x}{A \code{propsd} result.}
\item{digits}{Optional printing precision, otherwise following
\code{st.options}.}
\item{...}{Additional printing arguments.}
}
\details{
With both selected-sample arguments \code{NULL}, the selected sample equals
the complete reference. For a covariance matrix \eqn{B} describing
\eqn{n} sample occurrences, define \eqn{H=I-11'/n} and \eqn{C=HBH}.
Under a zero-mean Gaussian model the ordinary sample variance has exact
moments
\deqn{\mu_B=\mathrm{tr}(C)/(n-1),\qquad
v_B=2\sum_{ij}C_{ij}^2/(n-1)^2.}{
mu = trace(C)/(n-1); v = 2*sum(C^2)/(n-1)^2.
}
Moment matching gives Gamma shape \eqn{a=\mu_B^2/v_B} and scale
\eqn{b=v_B/\mu_B}. The expected standard deviation is approximated by
\deqn{e(B)=\sqrt{b}\,\Gamma(a+1/2)/\Gamma(a).}{
e(B) = sqrt(b)*Gamma(a+1/2)/Gamma(a).
}
For the selected multiset \eqn{s}, its covariance is the indexed matrix
\eqn{G_s=G[s,s]}, including duplicate rows and columns. The returned loss is
\deqn{1-e(G_s)/e(G).}
The sample-variance moments are exact under the model; the expected-SD formula
is generally an approximation. It is exact when the nonzero centered
covariance eigenvalues are equal. It is not the standard deviation of fitted
parental merits, group relatedness, or a mean of simulated SD ratios.

Multiplicities are evaluated without expanding a replicated matrix. For
\eqn{m=\sum_i k_i}, \eqn{w=k/m}, \eqn{r=Gw}, \eqn{q=w'Gw}, and
\eqn{D_{ij}=G_{ij}-r_i-r_j+q}, the selected moments are
\deqn{\mu_s=\frac{\sum_i k_iD_{ii}}{m-1},\qquad
v_s=\frac{2\sum_{ij}k_i k_jD_{ij}^2}{(m-1)^2}.}
An arbitrary scaling of all multiplicities can change the result because
ordinary sample variance uses \eqn{m-1}. The occurrence convention is therefore explicit.

A common positive covariance scale cancels from the loss. A common
random-intercept covariance cancels under centering. The implementation checks
symmetry and positive semidefiniteness, removes common intercept terms using
anchored contrasts, and scales covariance contrasts before evaluating moments
with extended-precision accumulators. Each selected active sample is anchored
and scaled directly from the original matrix; unused high-variance reference
parents therefore do not erase smaller selected contrasts. It uses log Gamma arithmetic; for shape
above 128 the log ratio correction uses
\eqn{-1/(8a)+1/(192a^3)-1/(640a^5)}. No replicated covariance workspace is made.
Reported moments use the original covariance units and can overflow or
underflow the R numeric range at extreme scales even when the loss and
expected SD remain well defined.

A zero-contrast reference is undefined and produces an error. A selected
sample with at least two identical parental occurrences has zero expected SD
and loss one. Negative losses are valid and are retained without clipping.
Input matrices, random-number states and thread settings are preserved.
}
\value{
A list of class \code{propsd}, containing:
\item{loss}{Computed proportional loss, possibly negative.}
\item{expected.sd}{Named vector \code{reference}, \code{selected} of
Gamma-approximated expected sample standard deviations in the square-root
units of the supplied covariance matrix.}
\item{n}{Named vector of reference and selected occurrence counts.}
\item{counts}{Selected multiplicities in complete reference order.}
\item{moments}{Data frame with rows \code{reference} and \code{selected};
columns \code{mean.variance}, \code{variance.variance}, \code{gamma.shape}
and \code{gamma.scale}. Gamma parameters are undefined and reported as
\code{NA} for a zero-contrast selected sample.}
\item{method}{Qualification identifying Gaussian sample-variance moments and
the two-moment Gamma approximation.}
}
\references{
Fern\'andez-Gonz\'alez, J., Metwally, S. M., and Isidro y S\'anchez, J.
(2026). MateR: a novel genomic mating framework.
\emph{Genetics} 232, iyag013. doi:10.1093/genetics/iyag013.
Equations 14--15 and Appendix C.
}
\seealso{\code{st.cross.plan}, \code{st.get.A}, \code{st.options}}
\examples{
data("v-tropmaize-vcf")
st.load.vcf.data(v.tropmaize.vcf, data.set = "propsd.example")
st.restrict.marker.data(NoAll.MAX = 2, MaMis.MAX = 0.05,
                        ExHet.MIN = 0.1, InMis.MAX = 0.1,
                        data.set = "propsd.example")
G <- st.get.A("propsd.example")/2
slots <- rownames(G)[c(1, 1, 2, 3)]
propsd(G, parents = slots)
counts <- tabulate(match(slots, rownames(G)), nbins = nrow(G))
propsd(G, counts = counts)
}
