\name{st.recode.ref}
\alias{st.recode.ref}
\title{Recode Marker Alleles Relative to a Reference Individual}
\description{
Recodes every marker genotype according to its allele sharing with one
reference individual.
}
\usage{
st.recode.ref(reference = 0, descending = FALSE, data.set = "default")
}
\arguments{
\item{reference}{Reference individual.  A character value is interpreted as an
individual identifier, whereas a positive whole-number numeric value is
interpreted as a one-based individual position in \code{data.set}.}
\item{descending}{Logical.  Selects the allele-sharing rule used for recoding;
see Details.}
\item{data.set}{Character string giving the SelectionTools marker data set to
recode.}
}
\details{
A character reference is looked up directly among the individual names in the
specified data set.  Thus \code{reference = "142"} means the individual named
\code{"142"}, whereas \code{reference = 142} means the individual at position
142.  The reference must identify exactly one existing individual; invalid
names and out-of-range numeric positions are rejected.

With \code{descending = FALSE} (the default), a genotype with the same two
alleles as the reference genotype is recoded to \code{1/1}; a genotype sharing
one reference allele is recoded to \code{1/2}; and a genotype sharing no
reference allele is recoded to \code{2/2}.  The reference individual is
therefore \code{1/1} wherever its source genotype is non-missing.

With \code{descending = TRUE}, a genotype is recoded to \code{1/1} when both
of its alleles can originate from the reference genotype, to \code{1/2} when
one allele can originate from the reference, and to \code{2/2} otherwise.
Missing source genotypes remain missing.

Whenever a genotype is recoded to \code{1/2}, the implementation stores the
reference-sharing state as \code{allele1 = 1} and \code{allele2 = 2}.
Consequently, this recoded \code{1/2} state is a summary of allele sharing
with the reference and must not be interpreted as inferred or preserved
gametic phase.
}
\value{
No useful value is intended.  The function modifies the selected marker data
set in place according to the reference-individual coding.
}
