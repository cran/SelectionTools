% copyright (C) 2001 Frisch/Maurer
%
\name{st.recode.hil}
\alias{st.recode.hil}
\title{st.recode.hil()}
\description{
  The function \emph{st.recode.hil} recodes the haplotype data for inbred populations.
}
\usage{
st.recode.hil(out.filename = "blocks", auxfiles = TRUE, data.set = "default")
}
\arguments{
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of data set that is used for recoding.}
}
\value{
Invisibly returns a \code{data.frame} with columns \code{Block}, \code{AlleleNr}, and \code{AlleleDef} describing recoded haplotype-block alleles when \code{auxfiles = TRUE} and the operation succeeds; otherwise \code{NULL}.
}

\details{The haplo block i of individual j at locus m is compared to the same haplo block
          of the other individuals. If the haplo blocks of different individuals are identical
          they get the same number. In extreme, there are as many different haplo blocks as individuals for
          a given locus. If a haplotype block would result in more than 254 distinct haplotype alleles,
          recoding stops with an error rather than merging or truncating haplotypes.\cr
          Using the function \emph{st.plot.ggt} after the function \emph{st.recode.hil} gives a good graphical
          overview, which haplo blocks are inherited from which parent. For plotting, the individuals are sorted
          by their haplotype frequency. Thus, the parents are usually plotted on top in yellow and blue (default values).
          The haplo blocks of the individuals are also plotted in yellow or blue, if they are completely identical
          to a parent's haplo block. If they are plotted in a different color, a cross-over did occur and thus,
          this haplo block is not identical one of the parent's haplo blocks.}

%\note{
%
%}
