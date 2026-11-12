\name{st.kernel.structure}
\alias{st.kernel.structure}
\title{Verify Kronecker structure in a relationship kernel}
\description{Detects a complete, ordered Kronecker product, or verifies a product
of two supplied parental relationship matrices on an explicit set of crosses.
Every matrix entry is checked. Detection does not replace a fitted kernel.}
\usage{
st.kernel.structure(K, factors = NULL, index = NULL,
                    dimensions = NULL, tolerance = 1e-10)
}
\arguments{
\item{K}{A finite, symmetric square numeric relationship matrix with nonnegative
diagonal. Row and column names, when supplied, must identify the same unique
levels; columns are aligned to row order. A failed structure search does not
certify positive semidefiniteness.}
\item{factors}{Optional list of exactly two positive semidefinite parental
relationship matrices. Their row and column names identify parental levels;
unnamed matrices receive positional labels \code{1,2,...}.}
\item{index}{When \code{factors} is supplied, a matrix or data frame with one row
per level of \code{K} and two columns in factor-list order. Entries are parental
names or one-based integer indices. Rows follow the row order of \code{K}.
Missing, out-of-range and unknown levels are rejected.}
\item{dimensions}{Optional pair of integer factor orders greater than one whose
product is the order of \code{K}. Without it, all nontrivial ordered splits are
tried until a verified factorization is found. Cannot be combined with
\code{factors} and \code{index}.}
\item{tolerance}{Maximum absolute reconstruction error divided by the maximum
absolute entry of \code{K}; between zero and \code{1e-6}. Default \code{1e-10}
allows numerical rounding. Approximate low-rank fitting is not performed.}
}
\details{
For a complete matrix in the corresponding order, \eqn{K=A\otimes B} has blocks
proportional to \eqn{B}. The function recovers candidate factors from principal
submatrices and verifies the full reconstruction without constructing a second
full Kronecker matrix. Unknown arbitrary permutations are not searched.

With supplied factors and indices it checks
\eqn{K_{ij}=c A_{d_i,d_j}B_{f_i,f_j}} for a positive scale \eqn{c}.
This also describes incomplete crosses and replicated crosses. A subset of a
Kronecker product need not itself admit a complete ordered factorization.
Zero kernels return no identified product because their factorization is
unidentified. Factors are not unique in general, including their reciprocal
scaling; the function reports one verified representation.

\code{st.lmm} automatically checks pairs of supplied relationship terms against
other relationship kernels, using compiled grouping indices. It requires a
consistent parental assignment for every target-kernel level. Rows excluded for
missing predictors do not establish mappings; missing phenotypes may do so.
Unmapped target levels are not guessed. It then tries complete ordered products
for unmatched kernels. The automatic mapped search checks at most 256 candidate
pairs per model; attributes \code{candidates.checked} and
\code{search.truncated} describe that search. The explicit function has no
candidate-pair limit.

Automatic results are available as \code{fit$model$kernel.products},
\code{fit$model$specification$kernel.products} and
\code{fit$engine$diagnostics$kernel.products}. The controls
\code{detect.kernel.products=FALSE} and \code{kernel.product.tol} disable or set
the tolerance of automatic detection. The metadata describes relationship
matrices, independently of estimated variance components or coefficient bases.
It is retained on serialization. In eligible scalar models, \code{st.lmm}
uses this result for a complete-crossing tensor solver or an observed-cell
spectral base with a smaller correction for other random terms. The selected
route retains the supplied covariance; it does not substitute a reconstructed
product. \code{kernel.solver=FALSE} keeps detection while disabling these
computational routes. Further details appear in \code{\link{st.lmm}} for eligibility, numerical
verification and diagnostic details.
}
\value{A list with \code{type} (\code{"none"}, \code{"kronecker"}, or
\code{"mapped.kronecker"}), \code{factors}, \code{index}, positive \code{scale},
\code{factor.dimensions}, \code{complete}, \code{max.absolute.error},
\code{relative.error}, \code{tolerance}, \code{checked.entries} and \code{reason}.
Indices refer to rows of the returned factors. \code{complete} requires every
parental pair exactly once. No-match results contain no factor or index guess.
Automatic detection returns a named list containing the verified matches only.}
\examples{
data("v-tropmaize-vcf", package = "SelectionTools")
marker.frame <- st.STvcf.to.dataframe(v.tropmaize.vcf)[1:12, c(1:6, 7:14)]
vcf <- st.dataframe.to.STvcf(marker.frame)
st.load.vcf.data(vcf, data.set = "maize.example")
relationship <- st.get.A("maize.example")
A <- relationship[1:2, 1:2]
B <- relationship[3:5, 3:5]
z <- st.kernel.structure(kronecker(A, B), dimensions = c(2, 3))
z$type
z$relative.error
}
\seealso{\code{\link{st.lmm}}}
