\name{st.get.A}
\alias{st.get.A}
\title{Constructs a genomic relationship matrix from a SelectionTools data set}
\description{Obtains the named marker design matrix from a SelectionTools marker data set and applies the existing VanRaden relationship
matrix constructor. Both operations are called from C.}
\usage{st.get.A(data.set = "default", coding = "0 1 2")}
\arguments{
\item{data.set}{Name of an existing SelectionTools marker data set, loaded
for example with \code{st.read.marker.data.df} or \code{st.load.vcf.data}.}
\item{coding}{Coding of the internal marker design matrix; passed unchanged
to \code{gs.make.A}. The default matches the standard internal 0/1/2 coding.}
}
\details{The compiled wrapper calls the existing named Z accessor, which builds
Z from the marker alleles if necessary, followed by the existing A constructor.
It does not implement a second marker converter or relationship formula.
For polymorphic biallelic input, the result is equivalent to
\code{gs.make.A(gs.get.Z(data.set, names=TRUE), coding=coding)}.
This relationship-matrix path also accepts fixed biallelic loci, which contribute zero after centering. An entirely fixed panel has zero relationship scale and is rejected. Multiallelic loci must be excluded. The separate legacy \code{gs.get.Z} interface retains its original marker requirements. The dedicated
\code{genomic_relationship_A} field caches the numeric relationship matrix.
The cache is owned by the compiled marker dataset in C memory and is independent
of returned R matrices. Removing or modifying an R result does not release or
alter this cache. compiled symbolic fits reuse it across calls. Its validity is checked
against marker content, individual and marker names and order, and coding;
performance-only changes do not invalidate it. Compaction clears the cache.
This exported accessor returns a copy, so editing the returned matrix cannot
change the internal cache. Fitted objects retain independent model snapshots;
removing the source dataset or clearing its cache does not invalidate a fit. No ridge or inverse is applied.}
\value{A numeric square relationship matrix. Individual identifiers label
both dimensions. Invalid inputs or construction failures produce an error.}
\examples{
data("v-tropmaize-vcf", package = "SelectionTools")
marker.frame <- st.STvcf.to.dataframe(v.tropmaize.vcf)[1:12, c(1:6, 7:14)]
vcf <- st.dataframe.to.STvcf(marker.frame)
st.load.vcf.data(vcf, data.set = "maize.example")
A <- st.get.A("maize.example")
A[1:4, 1:4]
}
\seealso{\code{\link{gs.get.Z}}, \code{\link{gs.make.A}},
\code{\link{st.read.marker.data.df}}, \code{\link{st.load.vcf.data}}}
