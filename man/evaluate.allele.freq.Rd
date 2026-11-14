\name{evaluate.allele.freq}
\alias{evaluate.allele.freq}
\title{Calculates or return allele frequencies for selected loci}

\description{
Calculates or return allele frequencies for selected loci.
}

\usage{
evaluate.allele.freq(PopName, Loci, missing=0)
}

\arguments{
   \item{PopName}{Name of the population.}
    \item{Loci}{A character string containing one locus name or
    multiple locus names separated by spaces. Character vectors
    of length greater than one are not supported.}
    
    \item{missing}{Integer controlling missing-allele handling.
    With \code{0}, missing alleles are included.
    With \code{1}, missing alleles are excluded.}
}
\value{
A \code{data.frame} containing the allele counts and relative frequencies.
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
evaluate.allele.freq("Maize", loci[1])
evaluate.allele.freq("Maize", paste(loci[1:2], collapse = " "))
}

\seealso{
Vector-based input of loci is supported by \code{\link{evaluate.allele.freq2}}.
}