\name{st.recode.ref}
\alias{st.recode.ref}
\title{st.recode.ref()}

\description{
  Compares the alleles of an individual with the alleles for a
  reference individual.
}

\usage{
st.recode.ref(reference = 0, descending = FALSE, data.set = "default")
}
\arguments{
   \item{reference}{Specifies which individual is the reference for recoding}
   \item{descending}{Check whether the alleles can originate from
   the reference indivdiual. Default: FALSE}
   \item{data.set}{Name of data set that is used for recoding}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of recoding marker data according to the reference-marker convention. Validation failures may return \code{NULL}.
}

\details{\code{descending=FALSE} (default):If an individuals has two alleles in
common with the reference, then it gets the genotype 1/1, for one allele
in common 1/2, and for non allele in common 2/2. The reference
individual is recoded to 1/1. Missing source genotypes remain missing.
\code{descending=TRUE}: If both alleles of an individual can originate
from the referende: 1/1, if one allele can originate from the reference:
1/2, if no allele can originate from the reference 2/2.
}
