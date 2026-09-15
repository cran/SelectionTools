% copyright (C) 2001 Frisch/Maurer
%
\name{st.recode.hbc}
\alias{st.recode.hbc}
\title{st.recode.hbc()}
\description{
  The function \emph{st.recode.hbc} recodes the haplotype data for backcross populations
  according to a given reference individual.
}
\usage{
st.recode.hbc(reference = 1, data.set = "default")
}
\arguments{
   \item{reference}{Specifies which individual is treated as reference for recoding. Usually the
                  recurrent parent is set to be the reference individual. Must be between 1 and
                  the number of individuals in the data set.}
   \item{data.set}{Name of data set that is used for recoding.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of recoding marker data according to the HBC reference convention. Validation failures may return \code{NULL}.
}

\details{If, at a given locus, an individual carries a haplotype that is also
         present in the reference genotype, then this haplotype is recoded to 1.
         If the haplotype is not present in the reference genotype, then it is
         recoded to 2.\cr
         Using the function \emph{st.plot.ggt} after the function \emph{st.recode.hbc} gives a good graphical
         overview, which haplo blocks are not inherited from the reference genotype (if reference
         individual = 1, the reference parent genome is plotted in yellow. The haplo
         blocks inherited from the donor parent are plotted in blue.}

%\note{
%
%}
