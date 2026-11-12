\name{st.hybrid.data.set}
\alias{st.hybrid.data.set}
\alias{st.hybrid.data.set.list}
\alias{print.st.hybrid.data.set.list}
\alias{st.hybrid.data.set.remove}
\alias{st.hybrid.data.set.remove.all}
\alias{st.set.hybrid.performance.data}
\alias{st.get.hybrid.performance.data}
\title{supported Hybrid Data Sets}
\description{Stores hybrid performance data with compiled references to two existing parental marker data sets, for symbolic GCA and SCA mixed models.}
\usage{
st.hybrid.data.set(markerdata1, markerdata2, hybrid.performance,
                   hybrid.data.set)
st.hybrid.data.set.list(digits = getOption("digits"))
\method{print}{st.hybrid.data.set.list}(x, digits = NULL, ...)
st.hybrid.data.set.remove(hybrid.data.set)
st.hybrid.data.set.remove.all()
st.set.hybrid.performance.data(hybrid.performance, hybrid.data.set)
st.get.hybrid.performance.data(hybrid.data.set)
}
\arguments{
\item{markerdata1,markerdata2}{Names of two existing, distinctly named marker data sets. Their markers are referenced, not copied.}
\item{hybrid.performance}{A nonempty data frame with parental ID columns named exactly like the two marker data sets. IDs must identify existing parental individuals. Response columns, missing responses, repeated observations, and additional model variables are retained.}
\item{hybrid.data.set}{One nonempty hybrid-data-set name. Marker and hybrid data sets share a namespace. Reusing an existing hybrid name replaces its performance and parent references after validation.}
\item{digits}{Number of decimal places for printing. A printing value of \code{NULL} uses \code{st.options()$print.digits}.}
\item{x}{A hybrid-data-set listing.}
\item{...}{Additional printing controls, including \code{top}, \code{bottom} and \code{char.per.line}; omitted values use \code{\link{st.options}}.}
}
\details{
The performance frame is stored as a snapshot. Changing the original R frame does not change the supported data. \code{st.set.hybrid.performance.data} replaces the stored frame and \code{st.get.hybrid.performance.data} retrieves a snapshot.

\code{st.lmm} accepts an existing response column, for example \code{fixed=GY~1}, with a quoted stored name such as \code{data="hybrids"}. With parental datasets named \code{flint} and \code{dent}, \code{random=~flint+dent} fits GCA only. Adding \code{flint:dent} fits SCA. \code{ucov=list(flint="A", dent="A", "flint:dent"="kro(a,a)")} specifies genomic GCA plus SCA. The product is one quoted string; \code{"kro(A,A)"} is also accepted. Prepared covariance lists work identically. Unquoted matrix names refer to actual R objects. No R function called \code{kro} or caller matrix called \code{A} is needed for the quoted symbolic requests.

Each marker data set owns a dedicated \code{genomic_relationship_A} cache. A hybrid data set owns \code{genomic_relationship_SCA}, a factorized representation containing the parental references and unique-cross mapping. A selected SCA covariance is cached only when requested. The full Cartesian Kronecker matrix is never constructed. Repeated observations of a cross share its SCA effect.

Marker cache validity is checked against marker content, individual and marker identifiers and order, and coding. Marker changes invalidate related SCA caches. A response-only update retains the relationship caches. Changing unique parental pairs invalidates the selected SCA cache. Compacting a parental data set clears disposable caches; they are reconstructed on demand. A referenced parental data set cannot be removed until its hybrid data sets have been removed. Removing a hybrid data set does not remove its parents. Existing fitted models retain their own snapshots.

Listings report parental names, observation and unique-cross counts, both parental cache flags, factorized and selected SCA cache flags, and estimated payload memory in MiB. Shared parental caches are counted only in their marker data sets. Factorized mappings and stored performance belong to hybrid \code{ram.core}; selected numeric SCA caches belong to hybrid \code{ram.work}.

The accessor \code{\link{st.lmm.components}(fit, "predictions")} returns the
stored hybrid prediction table with both parental identifiers.
For \code{prediction = "NAonly"}, \code{source.row} indexes the retained
performance data and identifies rows with missing responses. In a multivariate
fit, observed trait cells are left missing in this prediction table.
For \code{prediction = "Full"}, the table contains all combinations of the
two parental pools, with the second parent varying fastest and
\code{source.row = NA}. Full prediction requires non-parental fixed and
random covariates to be constant. The call
\code{predict(fit, newdata = ..., type = "conditional")} specifies other
covariate scenarios. Unrecorded SCA effects are predicted directly from the
parental factors and saved model scores, without inverting a genomic kernel.
The default mode is \code{"NAonly"}.
}
\value{
The constructor invisibly returns the supported name. The update function invisibly returns that name. The getter returns the stored performance frame as a copy. Listing returns a data frame of class \code{st.hybrid.data.set.list}. Removal functions invisibly return the number of hybrid data sets removed.
}
\seealso{\code{\link{st.lmm}}, \code{\link{st.get.A}}, \code{\link{st.data.set.list}}, \code{\link{st.compact.data.set}}}
