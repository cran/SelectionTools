\name{append.population}
\alias{append.population}
\title{append.population()}

\description{
Appends a copy of one population to another.
}

\usage{
append.population(NameP1, NameP2)
}

\arguments{
   \item{NameP1}{Character string naming the destination population.
   If it does not exist, it is created.}
   \item{NameP2}{Character string naming the source population.
   This population is not modified.}
}

\value{
Returns NameP1 and NameP2. Usually called for side effects.
}

\details{
The individuals of \code{NameP2} are copied to the end of
\code{NameP1}. Existing individuals in \code{NameP1} are retained,
and \code{NameP2} is not modified.

The source population must exist and contain at least one individual.
The source and destination population names must differ.
}

\examples{
data("v-tropmaize-vcf", package = "SelectionTools")
marker.frame <- st.STvcf.to.dataframe(v.tropmaize.vcf)[1:12, c(1:6, 7:14)]
vcf <- st.dataframe.to.STvcf(marker.frame)
st.load.vcf.data(vcf, data.set = "maize.example")
reset.all()
st.set.sim.gp(data.set = "maize.example")
st.set.sim.mp(data.set = "maize.example")
st.set.sim.pp("Maize", data.set = "maize.example")
loci <- vcf$ID
append.population("Combined", "Maize")
append.population("Combined", "Maize")
population.size.get("Combined")
}

\seealso{
\code{\link{population.append}}, which is called by this wrapper.
}